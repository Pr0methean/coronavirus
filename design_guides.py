from multiprocessing import Pool
from Bio import AlignIO, SeqIO
from datetime import datetime
from Bio.Seq import Seq
import itertools
import pickle
import shelve
import redis
import os

# length of the CRISPR guide RNA
K = 28
TARGETS_KEY = f"targets_{K}"
GOOD_TARGETS_KEY = f"good_targets_{K}"
# path to a fasta file with host sequences to avoid
HOST_FILE = "GCF_000001405.39_GRCh38.p13_rna.fna" # all RNA in human transcriptome
# HOST_FILE = "lung-tissue-gene-cds.fa" # just lungs
HOST_PATH = os.path.join("host", HOST_FILE)
# ending token for tries
END = "*"
# path to pickle / save the trie
REBUILD_TRIE = False
TRIE_PATH = "trie"
# path to alignment and id for the sequence to delete
TARGET_PATH = os.path.join("alignments", "HKU1+MERS+SARS+nCoV-Consensus.clu")
TARGET_ID = "nCoV"
# mismatch cutoff (e.g. how close must two kmers be to bind?)
CUTOFF = 2
# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 20
# paths to plasmid part fasta files
PROMOTER_PATH = os.path.join("parts", "pdpn_1_promoter.fa")
DR_SEQUENCE_PATH = os.path.join("parts", "crispr", "dr_sequence_bz_short.fa")
TAIL_PATH = os.path.join("parts", "tail.fa")
# path for output
OUTFILE_PATH = os.path.join("guides", "trie_guides.csv")

r = redis.Redis(host='localhost', port=6379)

# trie = shelve.open(TRIE_PATH)
trie = {}


# helpers
def all_equal(arr):
    return arr.count(arr[0]) == len(arr)


def read_fasta(fasta_path):
    record = SeqIO.read(handle=fasta_path, format="fasta")
    return record.seq.lower()


def write_fasta(fasta_path, sequences):
    SeqIO.write(sequences=sequences, handle=fasta_path, format="fasta")


def getKmers(sequence, k, step):
    for x in range(0, len(sequence) - k, step):
        yield sequence[x:x + k]


def index(kmer, haystack=trie):
    node = haystack
    for base in kmer:
        if base not in node:
            node[base] = {}
        node = node[base]
    node[END] = END


def _find(node, path, kmer, d):
    if not kmer:
        if len(path) is K:
            yield (path, d)
        return
    base, suffix = kmer[0], kmer[1:]
    for key in node:
        step = 1 if key is not base else 0
        if d + step > CUTOFF:
            return
        for result in _find(node[key], path + key, suffix, d + step):
            yield result


def find(kmer, haystack=trie):
    return _find(haystack, "", kmer, 0)


def host_has(kmer):
    matches = list(find(kmer))
    should_avoid = len(matches) > 0
    notice = "avoid" if should_avoid else "allow"
    print(notice, kmer, "matches", matches)
    return should_avoid


def make_hosts(input_path=HOST_PATH, db=r):
    if not REBUILD_TRIE:
        return
    trie.clear()
    with open(HOST_PATH, "r") as host_file:
        for rcount, record in enumerate(SeqIO.parse(host_file, "fasta")):
            for kmer in getKmers(record.seq.lower(), K, 1):
                kmer_string = str(kmer)
                r.sadd("hosts", kmer_string)
                index(kmer_string)
            print(rcount)


def make_targets(input_path=TARGET_PATH, target_id=TARGET_ID, db=r):

    alignment = AlignIO.read(input_path, "clustal")
    seq_ids = [seq.id for seq in alignment]
    index_of_target = seq_ids.index(target_id)
    alignment_length = alignment.get_alignment_length()
    conserved = [1 if all_equal(
        [seq[i] for seq in alignment]) else 0 for i in range(alignment_length)]
    for start in range(alignment_length - K):
        if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
            continue
        kmer = str(alignment[index_of_target][start:start + K].seq).lower()
        if "-" in kmer:
            continue
        n_conserved = sum(conserved[start:start + K])
        print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
        db.zadd(TARGETS_KEY, {kmer: n_conserved})
    most = db.zrevrangebyscore(TARGETS_KEY, 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer {most[0].decode()} has {int(most[1])} bases conserved in {seq_ids}")


def predict_side_effects(db=r):
    targets = db.zrevrangebyscore(TARGETS_KEY, 9001, 0)
    for target in targets:
        t = target.decode()
        should_avoid = host_has(t)
        if should_avoid:
            continue
        db.zadd(GOOD_TARGETS_KEY, {target: db.zscore(TARGETS_KEY, t)})
    with open(OUTFILE_PATH, "w+") as outfile:
        for k, good_target in enumerate(r.zrevrangebyscore(GOOD_TARGETS_KEY, 90, 0)):
            good_target_string = good_target.decode()
            print("good target", k, good_target_string)
            outfile.write(good_target_string + "\n")
    print(f"saved {db.zcard(GOOD_TARGETS_KEY)} good targets at {OUTFILE_PATH}")


if __name__ == "__main__":
    make_hosts()
    # test the trie lookup works
    for i in range(5):
        host_has(r.srandmember("hosts").decode())
    make_targets()
    predict_side_effects()
