# why?: predict side effects in biotech
# how?: identify K-length substrings of a target not present in a host
# what?: hamming distance cutoff
from cassandra.cluster import Cluster
from Bio import AlignIO, SeqIO
from itertools import product
import multiprocessing as mp
from tqdm import tqdm
import os


# cpus
CPUS = 12
# length of the CRISPR guide RNA
K = 28
# path to a fasta file with host sequences to avoid
HOST_FILE = "GCF_000001405.39_GRCh38.p13_rna.fna"  # human transcriptome
# HOST_FILE = "lung-tissue-gene-cds.fa" # lung protein codes
HOST_PATH = os.path.join("data", "host", HOST_FILE)
# path to pickle / save the index
REINDEX = True
# path to alignment and id for the sequence to delete
CORONAVIRUS_CONSENSUS = "HKU1+MERS+SARS+nCoV-Consensus.clu"
TARGET_PATH = os.path.join("data", "alignments", CORONAVIRUS_CONSENSUS)
TARGET_ID = "nCoV"
# mismatch cutoff (e.g. how close must two kmers be to bind?)
CUTOFF = 5
# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 21
# path for output
OUT_PATH = os.path.join("data", "guides", f"k{K}_cutoff{CUTOFF}_guides.csv")
# database
cluster = Cluster(contact_points=["172.28.1.1"])
global_session = cluster.connect()
# base -> options map
WILDCARD = {
    "n": ["a", "c", "t", "g"],
    "w": ["a", "t"],
    "a": ["a"],
    "c": ["c"],
    "t": ["t"],
    "g": ["g"]
}


def count_records(path):
    total = 0
    with open(path, "r") as file:
        for _ in SeqIO.parse(file, "fasta"):
            total = total + 1
    print(f"Total Records: {total}")
    return total


def count_kmers(path, k=K):
    total = 0
    for _ in generate_kmers(path, k):
        total = total + 1
    print(f"Total Kmers: {total}")
    return total


def get_kmers(record, k=K, stringify=1):
    rec = str(record.seq.lower()) if stringify else record
    if 'n' in rec or 'w' in rec:
        rec = [WILDCARD[base] for base in rec]
        for i in range(len(record) - k + 1):
            for option in product(*rec[i:i + k]):
                yield ''.join(option)
    else:
        for i in range(len(record) - k + 1):
            yield ''.join(rec[i:i+k])


def generate_kmers(path, k=K):
    with open(path, "r") as hostfile:
        for record in SeqIO.parse(hostfile, "fasta"):
            for kmer in get_kmers(record):
                yield kmer


SAVE_HOST = "insert into rna.hosts (kmer) values (?)"
SADD_UPDATE = "update rna.trie set next = next + ? where pre = ?"
NEXT = global_session.prepare("select * from rna.trie where pre = ?")


def init():
    global session, insert, update, contains, sadd, save, sismember
    session = cluster.connect()
    save = session.prepare(SAVE_HOST)
    update = session.prepare(SADD_UPDATE)


def _handle_kmer(kmer, k=K):
    session.execute_async(save, (kmer, ))
    for i in reversed(range(1, k)):
        session.execute_async(update, ({kmer[i]}, kmer[:i]))


def make_trie(path=HOST_PATH, cpus=CPUS, k=K):
    with mp.Pool(processes=CPUS, initializer=init) as pool:
        _ = list(
                tqdm(
                    pool.imap(_handle_kmer, generate_kmers(path, k=k)),
                    total=574817355)
            )


def _all_equal(arr):
    return arr.count(arr[0]) == len(arr)


ZADD = global_session.prepare(
    "insert into rna.targets (target, n, start, kmer, score, overlaps, host_has) values (?, ?, ?, ?, ?, ?, ?)")


def zadd(n, kmer, score, start, target='ncov', overlaps=True, host_has=True):
    global_session.execute_async(
        ZADD,
        (target, n, start, kmer, score, overlaps, host_has)
        )


def zrevrangebyscore(filter_overlaps=False):
    statement = "select n, kmer, score, start from rna.targets"
    if filter_overlaps:
        statement = statement + " where overlaps = False allow filtering"
    return global_session.execute(statement)


def overlap(s1, s2, k=K):
    if len(set(range(s1, s1+k)).intersection(range(s2, s2+k))) > 0:
        return True
    return False


def make_targets(path=TARGET_PATH, id=TARGET_ID, k=K,
                 offset_1=OFFSET_1, offset_2=OFFSET_2):
    alignment = AlignIO.read(path, "clustal")
    ids = [seq.id for seq in alignment]
    index_of_target = ids.index(id)
    alignment_length = alignment.get_alignment_length()
    conserved = [
        1 if _all_equal([seq[i] for seq in alignment]) else 0
        for i in range(alignment_length)
    ]
    n = 1
    for start in range(alignment_length - k + 1):
        kmer = str(alignment[index_of_target][start:start + k].seq).lower()
        offset_conserved = all(conserved[start + offset_1:start + offset_2])
        if "-" in kmer or not offset_conserved:
            continue
        n_conserved = sum(conserved[start:start + k])
        if n_conserved > 0:
            print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
            zadd(n, kmer, n_conserved, start)
            n = n + 1
    T = zrevrangebyscore()
    best = T.one()
    print(f"top {k}mer {best.kmer} has {int(best.score)} conserved in {ids}")
    no_overlaps = []
    for row in T:
        if any([overlap(other.start, row.start) for other in no_overlaps]):
            continue
        no_overlaps.append(row)
        zadd(row.n, row.kmer, row.score, row.start, overlaps=False)
    [print(r.kmer, r.score, r.start) for r in no_overlaps]
    print(f"{len(no_overlaps)} non-overlapping targets")
    no_overlaps_db = zrevrangebyscore(filter_overlaps=1)
    assert len(list(no_overlaps_db)) == len(no_overlaps)
    return no_overlaps


def _find(path, target, d, k=K, cutoff=CUTOFF):
    print("predict side effects", path, target, "distance", d)
    if not target:
        if len(path) is k:
            yield (path, d)
        return
    target_base, suffix = target[0], target[1:]
    pre = path if path else "root"
    for host_base in global_session.execute(NEXT, (pre,)).one().next:
        step = 1 if host_base != target_base else 0
        if d + step >= CUTOFF:
            return
        for result in _find(path + host_base, suffix, d + step, k, cutoff):
            yield result


def _host_has(target, cutoff=CUTOFF, k=K):
    matches = list(_find("", target, 0, cutoff=cutoff, k=k))
    print(target, "matches", matches)
    return len(matches) > 0


def predict_side_effects(k=K, cutoff=CUTOFF):
    for t in zrevrangebyscore(filter_overlaps=1):
        if _host_has(t.kmer, cutoff=cutoff):
            continue
        print(
            f"good target {t.kmer} is {int(t.score)}/{k} conserved ("
            f"{int((t.score / k) * 100)}%)")
        zadd(t.n, t.kmer, t.score, t.start, overlaps=False, host_has=False)


if __name__ == "__main__":
    make_trie()
    make_targets()
    predict_side_effects()
