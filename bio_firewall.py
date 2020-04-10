# why?: predict side effects in biotech
# how?: identify K-length substrings of a target not present in a host
# what?: hamming distance cutoff
from cassandra.cluster import Cluster
import multiprocessing as mp
from itertools import product
import os

import frozendict as frozendict
from Bio import AlignIO, SeqIO
from redis import Redis
from tqdm import tqdm

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
CUTOFF = 4
# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 21
# path for output
OUT_PATH = os.path.join("data", "guides", f"k{K}_cutoff{CUTOFF}_guides.csv")
# cache
REDIS_ARGS = frozendict.frozendict()
r = Redis(**REDIS_ARGS)
# database
cluster = Cluster(contact_points=["127.0.1.1"])
session, sadd_stmt, sismember_stmt = None, None, None
# base -> options map
WILDCARD = {
    "n": ["a", "c", "t", "g"],
    "w": ["a", "t"],
    "a": ["a"],
    "c": ["c"],
    "t": ["t"],
    "g": ["g"]
}
SADD_INSERT = "insert into rna.trie (pre) values (?) if not exists"
SADD_UPDATE = "update rna.trie set next = next + ? where pre = ?"
SISMEMBER = "select * from rna.trie where pre = ? and next contains ? ALLOW FILTERING"


def init():
    global session, sadd, sadd_stmt, sismember, sismember_stmt
    session = cluster.connect()
    sadd_insert = session.prepare(SADD_INSERT)
    sadd_update = session.prepare(SADD_UPDATE)
    sismember_stmt = session.prepare(SISMEMBER)

    def sadd(k, v):
        insert_result = session.execute(sadd_insert, [k])
        update_result = session.execute(sadd_update, [{v}, k])
        return insert_result and update_result

    def sismember(k, v):
        return session.execute(sismember_stmt, [k, v]).one() is not None


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


def _handle_kmer(kmer, k=K):
    if sismember("hosts", kmer):
        return
    sadd("hosts", kmer)
    for i in reversed(range(1, k)):
        prefix, base = kmer[:i], kmer[i]
        sadd(prefix, base)
        if sismember(prefix[:-1], prefix[-1]):
            break
    return


def make_hosts(path=HOST_PATH, cpus=CPUS, k=K, redis_args=REDIS_ARGS,
               reindex=REINDEX):
    if reindex:
        total = 574817355
        with mp.Pool(processes=CPUS - 2, initializer=init) as pool:
            _ = list(
                    tqdm(
                        pool.imap(_handle_kmer, generate_kmers(path, k=k)),
                        total=total)
                )
    r = Redis()
    return [x.decode() for x in r.smembers("hosts")]


def _find(path, target, d, k=K, db=r, cutoff=CUTOFF):
    if not target:
        if len(path) is k:
            yield (path, d)
        return
    target_base, suffix = target[0], target[1:]
    node = path if path else "root"
    for option in db.smembers(node):
        host_base = option.decode()
        step = 1 if host_base != target_base else 0
        if d + step > CUTOFF:
            return
        for result in _find(path + host_base, suffix, d + step, k, db, cutoff):
            yield result


def _host_has(target, cutoff=CUTOFF, k=K):
    matches = list(_find("", target, 0, cutoff=cutoff, k=k))
    print(target, "matches", matches)
    return len(matches) > 0


def _all_equal(arr):
    return arr.count(arr[0]) == len(arr)


def make_targets(path=TARGET_PATH, id=TARGET_ID, k=K, db=r, offset_1=OFFSET_1,
                 offset_2=OFFSET_2):
    targets_key = f"targets_{k}"
    alignment = AlignIO.read(path, "clustal")
    ids = [seq.id for seq in alignment]
    index_of_target = ids.index(id)
    alignment_length = alignment.get_alignment_length()
    conserved = [
        1 if _all_equal([seq[i] for seq in alignment]) else 0
        for i in range(alignment_length)
    ]
    for start in range(alignment_length - k + 1):
        kmer = str(alignment[index_of_target][start:start + k].seq).lower()
        if "-" in kmer or not all(conserved[start + offset_1:start + offset_2]):
            continue
        n_conserved = sum(conserved[start:start + k])
        if n_conserved > 0:
            print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
            db.zadd(targets_key, {kmer: n_conserved})
    targets = db.zrevrangebyscore(targets_key, 9001, 0, withscores=True)
    T = [(t[0].decode(), t[1]) for t in targets]
    print(f"top {k}mer {T[0][0]} has {T[0][1]} bases conserved in {ids}")
    return T


def predict_side_effects(k=K, cutoff=CUTOFF, db=r):
    if not db.exists("hosts") and db.exists("targets"):
        raise Exception("hosts and targets required to predict side effects")
    targets_key, good_targets_key = f"targets_{k}", f"good_targets_{k}"
    targets = db.zrevrangebyscore(targets_key, 9001, 0)
    for target in targets:
        t = target.decode()
        if _host_has(t, cutoff=cutoff):
            continue
        db.zadd(good_targets_key, {t: db.zscore(targets_key, t)})
    good_targets = db.zrevrangebyscore(
        good_targets_key, 90, 0, withscores=True)
    GT = [(gt[0].decode(), gt[1]) for gt in good_targets]
    [print(
        f"good target {gt[0]} has {int(gt[1])} of {k} bases conserved ("
        f"{int((gt[1] / k) * 100)}%)")
        for gt in GT]
    with open(OUT_PATH, "w+") as outfile:
        outfile.writelines([f"{gt[0]}, {gt[1]}\n" for gt in GT])
    return GT


if __name__ == "__main__":
    make_hosts()
    # make_targets()
    # predict_side_effects()
