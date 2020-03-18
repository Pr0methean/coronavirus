from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import itertools
import redis
import os
from datetime import datetime

# length of the CRISPR guide RNA
K = 28
# path to a fasta file with host sequences to avoid
HOST_PATH = os.path.join("host", "lung-tissue-gene-cds.fa")
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
OUTFILE_PATH = os.path.join(
    "guides", "SARS-nCoV-2_consensus_conserved_watson_crick_guides_RNA.csv")

r = redis.Redis(host='localhost', port=6379)
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


# TODO: rewrite with "index" function
def make_hosts(input_path=HOST_PATH, db=r):
    rcount, kcount = 0, 0
    with open(input_path, "r") as host_file:
        for record in SeqIO.parse(host_file, "fasta"):
            rcount = rcount + 1
            for kmer in getKmers(record.seq.lower(), K, 1):
                kcount = kcount + 1
                print(rcount, kcount, kmer)
                db.sadd("hosts", str(kmer))


def make_targets(input_path=TARGET_PATH, target_id=TARGET_ID, db=r):
    alignment = AlignIO.read(input_path, "clustal")
    sequence_ids = [seq.id for seq in alignment]
    index_of_target = sequence_ids.index(target_id)
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
        db.zadd("targets", {kmer: n_conserved})
    most_conserved_kmer = db.zrevrangebyscore(
        "targets", 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer is {most_conserved_kmer[0].decode()} with {int(most_conserved_kmer[1])} bases conserved between {sequence_ids}")


# TODO: rewrite with "host_has" function
def predict_side_effects(db=r):
    targets = db.zrevrangebyscore("targets", 9001, 0)
    hosts = db.smembers("hosts")
    for target in targets:
        for host in hosts:
            d = sum([0 if target[n] is host[n] else 1 for n in range(K)])
            print(f"d({target.decode()}, {host.decode()}) = {d}")
            if d < CUTOFF:
                print("found potential mismatch:",
                      target.decode(), host.decode())
                break
        print("no side effects found for: ", target)
        db.zadd("good_targets", {target: db.zscore("targets", target)})


def make_plasmids(prefix_path=PROMOTER_PATH, suffix_path=DR_SEQUENCE_PATH, db=r):
    prefix = read_fasta(prefix_path)
    suffix = read_fasta(suffix_path)
    # tail = read_fasta(TAIL_PATH)
    good_targets = db.zrevrangebyscore("good_targets", 9001, 0)
    timestamp = datetime.now()
    for i, target in enumerate(good_targets):
        guide = target.reverse_complement()
        score = db.zscore(target)
        title = f"{i}_{guide}_{score}_{timestamp}"
        transcripts = [prefix, guide, suffix]
        plasmid = "".join(transcripts)
        write_fasta(title, plasmid)


if __name__ == "__main__":
    # make_hosts()
    # make_targets()
    # predict_side_effects()
    make_plasmids()
