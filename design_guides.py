from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import itertools
import redis
import os
from datetime import datetime
from multiprocessing import Pool

r = redis.Redis(host='localhost', port=6379)

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
OUTFILE_PATH = os.path.join("guides", "SARS-nCoV-2_consensus_conserved_watson_crick_guides_RNA.csv")


# helper functions
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


def handle_host_record(record):
    for kmer in getKmers(record.seq.lower(), K, 1):
        r.sadd("hosts", str(kmer))


def make_hosts():
    rcount, kcount = 0, 0
    with open(HOST_PATH, "r") as host_file:
        for record in SeqIO.parse(host_file, "fasta"):
            rcount = rcount + 1
            for kmer in getKmers(record.seq.lower(), K, 1):
                kcount = kcount + 1
                print(rcount, kcount, kmer)
                r.sadd("hosts", str(kmer))


def make_targets():
    alignment = AlignIO.read(TARGET_PATH, "clustal")
    sequence_ids = [seq.id for seq in alignment]
    index_of_target = sequence_ids.index(TARGET_ID)
    alignment_length = alignment.get_alignment_length()
    conserved = [1 if all_equal([seq[i] for seq in alignment]) else 0 for i in range(alignment_length)]
    for start in range(alignment_length - K):
        if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
            continue
        kmer = str(alignment[index_of_target][start:start + K].seq).lower()
        if "-" in kmer:
            continue
        n_conserved = sum(conserved[start:start + K])
        print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
        r.zadd("targets", {kmer: n_conserved})
        r.set(f"{kmer}:start", start)
    most_conserved_kmer = r.zrevrangebyscore("targets", 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer is {most_conserved_kmer[0].decode()} with {int(most_conserved_kmer[1])} bases conserved between {sequence_ids}")


def predict_side_effects():
    targets = r.zrevrangebyscore("targets", 9001, 0)
    hosts = r.smembers("hosts")
    for target in targets:
        for host in hosts:
            d = sum([0 if target[n] is host[n] else 1 for n in range(K)])
            print(f"d({target.decode()}, {host.decode()}) = {d}")
            if d < CUTOFF:
                print("found potential mismatch:", target.decode(), host.decode())
                break
        print("no side effects found for: ", target)
        r.zadd("good_targets", {target: r.zscore("targets", target)})


# def save_results():
#    with open(OUTFILE_PATH, "w+") as outfile:
#        outfile.write("\n".join(list(guides)))


def make_plasmids():
    pol3_promoter = read_fasta(PROMOTER_PATH)
    dr_sequence = read_fasta(DR_SEQUENCE_PATH)
    # tail = read_fasta(TAIL_PATH)
    good_targets = r.zrevrangebyscore("good_targets", 9001, 0)
    timestamp = datetime.now()
    for i, target in enumerate(good_targets):
        guide = target.reverse_complement()
        score = r.zscore(target)
        title = f"{i}_{guide}_{score}_{timestamp}"
        transcripts = [pol3_promoter, guide, dr_sequence]
        plasmid = "".join(transcripts)
        write_fasta(title, plasmid)


if __name__ == "__main__":
    # make_hosts()
    # make_targets()
    # predict_side_effects()
    make_plasmids()

