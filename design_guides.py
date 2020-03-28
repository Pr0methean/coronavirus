import os

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from redis import Redis
from multiprocessing import Pool


def bytesu(string):
    return bytes(string, "UTF-8")


# length of the CRISPR guide RNA
K = 28
# path to a fasta file with host sequences to avoid
HOST_FILE = "GCF_000001405.39_GRCh38.p13_rna.fna"  # all RNA in human transcriptome
# HOST_FILE = "lung-tissue-gene-cds.fa" # just lungs
HOST_PATH = os.path.join("host", HOST_FILE)
# ending token for tries
END = bytesu("*")
EMPTY = bytesu('')
# path to pickle / save the trie
REBUILD_TRIE = True
TRIE_PATH = "trie"
# path to alignment and id for the sequence to delete
TARGET_PATH = os.path.join("alignments", "HKU1+MERS+SARS+nCoV-Consensus.clu")
TARGET_ID = "nCoV"
# mismatch cutoff (e.g. how close must two kmers be to bind?)
CUTOFF = 2
# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 21
# paths to plasmid part fasta files
PROMOTER_PATH = os.path.join("parts", "pdpn_1_promoter.fa")
DR_SEQUENCE_PATH = os.path.join("parts", "crispr", "dr_sequence_bz_short.fa")
TAIL_PATH = os.path.join("parts", "tail.fa")
# path for output
OUTFILE_PATH = os.path.join("guides", "trie_guides.csv")


global my_db
global my_tdb
global my_k


# helpers
def all_equal(arr):
    return arr.count(arr[0]) == len(arr)


def read_fasta(fasta_path: str) -> Seq:
    record = SeqIO.read(handle=fasta_path, format="fasta")
    return record.seq.lower()


def write_fasta(fasta_path: str, sequences):
    SeqIO.write(sequences=sequences, handle=fasta_path, format="fasta")


def getKmers(sequence: str, k: int, step: int):
    for x in range(0, len(sequence) - k + 1, step):
        yield sequence[x:x + k]


def index(kmer: str, db: Redis):
    if db.get(kmer):
        return
    db.sadd(kmer, '*')
    for x in reversed(range(0, len(kmer))):
        prefix = kmer[:x]
        if not db.sadd(prefix, bytesu(kmer[x])):
            return


def _find(path: bytes, kmer: bytes, d: int, db: Redis, max_mismatches=CUTOFF, k=K):
    if not kmer:
        if len(path) is k:
            yield (path, d)
        return
    branches = db.sscan_iter(str(path, "UTF-8"))
    if not branches:
        return
    base, suffix = kmer[0], kmer[1:]
    for key in branches:
        step = 1 if key is not base else 0
        if d + step > max_mismatches:
            return
        for result in _find(path + key, suffix, d + step, db, max_mismatches, k):
            yield result


def find(kmer: str, tdb: Redis, max_mismatches=CUTOFF, k=K):
    return _find(b'', bytesu(kmer), 0, tdb, max_mismatches, k)


def host_has(kmer: str, tdb: Redis, max_mismatches=CUTOFF, k=K):
    matches = list(find(kmer, tdb, max_mismatches, k))
    should_avoid = len(matches) > 0
    notice = "avoid" if should_avoid else "allow"
    print(notice, kmer, "matches", matches)
    return should_avoid


def index_record(record_pair):
    (record_count, record) = record_pair
    for kmer in getKmers(record.seq.lower(), my_k, 1):
        kmer_string = str(kmer)
        my_db.sadd("hosts", bytesu(kmer_string))
        index(kmer_string, my_tdb)
    print(record_count)


def init_subprocess(db, tdb, k):
    global my_db, my_tdb, my_k
    my_db = db()
    my_tdb = tdb()
    my_k = k


def make_hosts(input_path: str, db, tdb, k=K):
    with Pool(initializer=lambda: init_subprocess(db, tdb, k)) as pool, open(input_path, "r") as host_file:
        records = enumerate(SeqIO.parse(host_file, "fasta"))
        pool.map(index_record, records)


def make_targets(db, target_path=TARGET_PATH, target_id=TARGET_ID, k=K):
    targets_key = f"targets_{k}"
    alignment = AlignIO.read(target_path, "clustal")
    seq_ids = [seq.id for seq in alignment]
    index_of_target = seq_ids.index(target_id)
    alignment_length = alignment.get_alignment_length()
    conserved = conserved_in_alignment(alignment, alignment_length)
    for start in range(alignment_length - k + 1):
        kmer, n_conserved = count_conserved(alignment, conserved, index_of_target, start, k)
        if n_conserved > 0:
            print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
            db.zadd(targets_key, {kmer: n_conserved})
    most = db.zrevrangebyscore(targets_key, 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer {most[0].decode()} has {int(most[1])} bases conserved in {seq_ids}")


def conserved_in_alignment(alignment, alignment_length):
    return [1 if all_equal(
        [seq[i] for seq in alignment]) else 0 for i in range(alignment_length)]


def count_conserved(alignment, conserved, index_of_target, start, k=K):
    if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
        return "", 0
    else:
        kmer = str(alignment[index_of_target][start:start + k].seq).lower()
        if "-" in kmer:
            return "", 0
        n_conserved = sum(conserved[start:start + k])
    return kmer, n_conserved


def predict_side_effects(db: Redis, out_path: str, tdb: Redis, k=K, max_mismatches=CUTOFF):
    targets_key = f"targets_{k}"
    good_targets_key = f"good_targets_{k}"
    targets = db.zrevrangebyscore(targets_key, 9001, 0)
    for target in targets:
        t = target.decode()
        should_avoid = host_has(t, tdb, max_mismatches, k)
        if should_avoid:
            continue
        db.zadd(good_targets_key, {t: db.zscore(targets_key, t)})
    with open(out_path, "w+") as outfile:
        for k, good_target in enumerate(db.zrevrangebyscore(good_targets_key, 90, 0)):
            good_target_string = good_target.decode()
            print("good target", k, good_target_string)
            outfile.write(good_target_string + "\n")
    print(f"saved {db.zcard(good_targets_key)} good targets at {out_path}")


if __name__ == "__main__":
    get_r = lambda: Redis(host='localhost', port=6379, db=0)
    r = get_r()
    get_r_for_trie = lambda: Redis(host='localhost', port=6379, db=1)
    r_for_trie = get_r_for_trie()
    if REBUILD_TRIE:
        make_hosts(HOST_PATH, db=get_r, tdb=get_r_for_trie)
    # test the trie lookup works
    for i in range(5):
        host_has(r.srandmember("hosts").decode(), tdb=r_for_trie)
    make_targets(db=r)
    predict_side_effects(db=r, tdb=r_for_trie, out_path=OUTFILE_PATH)
