import os

import plyvel
import redis
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq


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

r = None
leveldb = None


# helpers
def all_equal(arr):
    return arr.count(arr[0]) == len(arr)


def read_fasta(fasta_path: str) -> Seq:
    record = SeqIO.read(handle=fasta_path, format="fasta")
    return record.seq.lower()


def write_fasta(fasta_path: str, sequences):
    SeqIO.write(sequences=sequences, handle=fasta_path, format="fasta")


def getKmers(sequence: str, k: int, step: int):
    for x in range(0, len(sequence) - k, step):
        yield sequence[x:x + k]


def index(kmer: str, db=leveldb):
    kmer_bytes = bytesu(kmer + '*')
    with db.write_batch() as wb:
        for x in reversed(range(1, len(kmer_bytes))):
            prefix = kmer_bytes[:x]
            old_value = db.get(prefix, EMPTY)
            if kmer_bytes[x] in old_value:
                return  # this prefix is shared with a kmer that's already indexed
            else:
                new_value = bytes(sorted([*old_value, kmer_bytes[x]]))
                wb.put(prefix, new_value)


def _find(path, kmer, d, db, max_mismatches):
    if not kmer:
        if len(path) is K:
            yield (path, d)
        return
    base, suffix = kmer[0], kmer[1:]
    for key in db.get(bytesu(path), bytes()):
        step = 1 if key is not bytesu(base) else 0
        if d + step > max_mismatches:
            return
        for result in _find(path + key, suffix, d + step, db, max_mismatches):
            yield result


def find(kmer, db=leveldb, max_mismatches=CUTOFF):
    return _find("", kmer, 0, db, max_mismatches)


def host_has(kmer, db=leveldb):
    matches = list(find(kmer, db))
    should_avoid = len(matches) > 0
    notice = "avoid" if should_avoid else "allow"
    print(notice, kmer, "matches", matches)
    return should_avoid


def make_hosts(input_path=HOST_PATH, db=r, ldb=leveldb):
    if not REBUILD_TRIE:
        return
    with open(input_path, "r") as host_file:
        for rcount, record in enumerate(SeqIO.parse(host_file, "fasta")):
            for kmer in getKmers(record.seq.lower(), K, 1):
                kmer_string = str(kmer)
                db.sadd("hosts", kmer_string)
                index(kmer_string, ldb)
            print(rcount)


def make_targets(db=r, target_path=TARGET_PATH, target_id=TARGET_ID):
    alignment = AlignIO.read(target_path, "clustal")
    seq_ids = [seq.id for seq in alignment]
    index_of_target = seq_ids.index(target_id)
    alignment_length = alignment.get_alignment_length()
    conserved = conserved_in_alignment(alignment, alignment_length)
    for start in range(alignment_length - K):
        kmer, n_conserved = count_conserved(alignment, conserved, index_of_target, start)
        if n_conserved > 0:
            print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
            db.zadd("targets", {kmer: n_conserved})
    most = db.zrevrangebyscore("targets", 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer {most[0].decode()} has {int(most[1])} bases conserved in {seq_ids}")


def conserved_in_alignment(alignment, alignment_length):
    return [1 if all_equal(
        [seq[i] for seq in alignment]) else 0 for i in range(alignment_length)]


def count_conserved(alignment, conserved, index_of_target, start):
    if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
        kmer = ""
        n_conserved = 0
    else:
        kmer = str(alignment[index_of_target][start:start + K].seq).lower()
        if "-" in kmer:
            n_conserved = 0
        else:
            n_conserved = sum(conserved[start:start + K])
    return kmer, n_conserved


def predict_side_effects(db=r, out_path=OUTFILE_PATH):
    targets = db.zrevrangebyscore("targets", 9001, 0)
    for target in targets:
        t = target.decode()
        should_avoid = host_has(t)
        if should_avoid:
            continue
        db.zadd("good_targets", {target: db.zscore("targets", t)})
    with open(out_path, "w+") as outfile:
        for k, good_target in enumerate(db.zrevrangebyscore("good_targets", 90, 0)):
            good_target_string = good_target.decode()
            print("good target", k, good_target_string)
            outfile.write(good_target_string + "\n")
    print(f"saved {db.zcard('good_targets')} good targets at {out_path}")


if __name__ == "__main__":
    r = redis.Redis(host='localhost', port=6379)
    leveldb = plyvel.DB("db/", create_if_missing=True)
    make_hosts()
    # test the trie lookup works
    for i in range(5):
        host_has(r.srandmember("hosts").decode())
    make_targets()
    predict_side_effects()
    leveldb.close()
