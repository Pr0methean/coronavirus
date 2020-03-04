# conda create --name bio python=3.7.1 biopython redis
# conda activate bio

from Bio.Seq import Seq
from Bio import SeqIO
import itertools
import redis
import os

r = redis.Redis(host='localhost', port=6379)

# length of the guide RNA for CRISPR Cas13
K = 30

# selected non-watson-crick base pair map (cDNA right now)
WOBBLES = {
  "a": ["t"],
  "c": ["g"],
  "g": ["t", "c"],
  "t": ["g", "a"],
  "y": ["g", "a"], # y is c or t, c binds g and t binds g or a (in RNA form)
  "n": ["a", "c", "t", "g"]
}

# redis key for the set of whitelist (host/off-target) kmers
whitelist_key = "whitelist:lung:cds"

# sequence fasta for the target
target_fasta_file_name = "4-march-2020_57xSARS-nCoV-2_consensus.fa"
FASTA_PATH = os.path.join("blacklist", target_fasta_file_name)

# path for output
OUTFILE_PATH = os.path.join("guides", "SARS-nCoV-2_consensus_watson_crick_guides.csv")

# helper functions
def read_fasta(fasta_path):
  record = SeqIO.read(fasta_path, "fasta")
  return record.seq.lower()

def getKmers(sequence, k, step):    
  for x in range(0, len(sequence) - k, step):
    yield sequence[x:x+k]

def wobble(sequence, label1, label2):
  for seq in itertools.product(*[WOBBLES[char] for char in sequence]):
    string = "".join(seq)
    print(label1, ":", sequence, "-->", label2, ":", string)
    yield string

def design_guides():
  ncov_consensus = read_fasta(FASTA_PATH)
  kmers = getKmers(str(ncov_consensus), 30, 1)
  guides, count, hits = set(), 0, 0
  for kmer in kmers:
    count = count + 1
    print(count, kmer)
    if r.sismember(whitelist_key, kmer):
      hits = hits + 1
      continue
    guides.add(str(Seq(kmer).reverse_complement()))
  print("guides", guides, "count", count, "hits", hits, "percent:", hits / count * 100)
  with open(OUTFILE_PATH, "w+") as outfile:
    outfile.write("\n".join(list(guides)))

def design_guides_with_wobble():
  # read the coronavirus consensus
  ncov_consensus = read_fasta(FASTA_PATH)
  # generate kmers
  kmers = getKmers(str(ncov_consensus), 30, 1)
  # create guides for each kmer ... and targets for each guide
  # if a guide has any targets in the whitelist, we skip it
  i = 0
  guides = set()
  for kmer in kmers:
    print(i, kmer)
    i = i + 1
    if r.sismember(whitelist_key, kmer):
      continue
    for guide in wobble(kmer, "kmer", "guide"):
      try:
        targets = set() 
        for target in wobble(guide, "guide", "target"):
          if r.sismember(whitelist_key, target):
            print("HIT!")
            raise Exception("guide target in whitelist")
          targets.add(target)
        guides.add(guide)
        print("guide:", guide, "targets:", targets)
      except Exception as e:
        print(e)
        pass
  print("guides:", guides)
  with open(os.path.join("guides", target_fasta_file_name), "w") as outfile:
    outfile.write("\n".join(list(guides)))  

if __name__ == "__main__":
  design_guides()

# def find_wobbles(sequence, nucleotides=['g', 't', 'y', 'n']):
#   return [(i, sequence[i]) for i in range(len(sequence)) if sequence[i] in nucleotides]

# build a whitelist
# count = 0
# with open(os.path.join("whitelist", "lung-tissue-gene-cds.fa"), "r") as human_transcriptome:
#   for record in SeqIO.parse(human_transcriptome, "fasta"):
#     print(record.id)
#     kmers = getKmers(record.seq.lower(), K, 1)
#     for kmer in kmers:
#       # count = count + 1
#       # print(count, kmer)
#       r.sadd("whitelist:lung:cds", str(kmer))

# save the whitelist for later
# with open(os.path.join("whitelist", f"lung-cds-{K}mers.txt"), "w") as outfile:
#   outfile.write("\n".join(r.smembers('whitelist:lung:cds')))

# conserved substrings found manually in a Clustal alignment w/ ctrl-f ninjitsu
# CONSERVED = [
#   "ccttataattcacagaa",
#   "tgtggtagtgttggtt",
#   "tatgggttgggatta",
#   "tttaaatattggga",
#   "ggtttctttaagga",
#   "actcaaatgaatct",
#   "acggtttcgtccg",
#   "cttttcaaactgt",
#   "attcacagacttc",
#   "aattggaattgt",
#   "cgtgtacgccaa",
#   "ttcttgttaaca",
#   "gtggtcattcaa"
# ]

# def checkConserved(kmer):
#   conserved = any(fragment in kmer for fragment in CONSERVED)
#   return conserved

# with open(os.path.join("blacklist", "sars-cov-2.fasta"), "r") as coronavirus_genome:
#   for record in SeqIO.parse(coronavirus_genome, "fasta"):
#     covid_kmers = getKmers(record.seq.lower(), K, 1)
#     for kmer in covid_kmers:
#       if checkConserved(kmer):
#         print(kmer.transcribe().reverse_complement())

# def get_max_uppercase_run_from_string(s):
#     # construct a list of all the uppercase segments in your string
#     list_of_uppercase_runs = re.findall(r"[A-Z]+", s)
#     # sort by length
#     sorted_by_length = sorted(list_of_uppercase_runs, key=len, reverse=False)
#     [print(len(uppercase_run), uppercase_run) for uppercase_run in sorted_by_length]
#     # find out what the longest string is in your list
#     longest_string = sorted_by_length[-1]
#     # return the length of this string to the user
#     return longest_string

# with open("/Users/bionhoward/hax/coronavirus/blacklist/meta-consensus-nCoV-consensus+MERS+HKU1+SARS.fa", "r") as consensus:
#   for record in SeqIO.parse(consensus, "fasta"):
#     get_max_uppercase_run_from_string(str(record.seq))