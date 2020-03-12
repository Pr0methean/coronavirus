from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import itertools
import redis
import os

r = redis.Redis(host='localhost', port=6379)

# length of the guide RNA for CRISPR Cas13
K = 28

# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 20

# path to the alignment file to find conserved kmers
ALIGNMENT_PATH = os.path.join("alignments", "HKU1+MERS+SARS+nCoV-Consensus.clu")

# fasta id for the specific subsequence of the alignment to target
TARGET_ID = "nCoV"

# redis key for the set of whitelist (host/off-target) kmers
whitelist_key = "whitelist:lung:cds"

# path to a fasta file with host sequences to avoid
whitelist_file_path = os.path.join("whitelist", "lung-tissue-gene-cds.fa")

# sequence fasta for the target
target_fasta_file_name = "4-march-2020_57xSARS-nCoV-2_consensus.fa"
FASTA_PATH = os.path.join("blacklist", target_fasta_file_name)

# path for output
OUTFILE_PATH = os.path.join("guides", "SARS-nCoV-2_consensus_conserved_watson_crick_guides_RNA.csv")

# conserved substrings found manually in a Clustal alignment w/ ctrl-f ninjitsu
CONSERVED = [
  "ccttataattcacagaa",
  "tgtggtagtgttggtt",
  "tatgggttgggatta",
  "tttaaatattggga",
  "ggtttctttaagga",
  "actcaaatgaatct",
  "acggtttcgtccg",
  "cttttcaaactgt",
  "attcacagacttc",
  "aattggaattgt",
  "cgtgtacgccaa",
  "ttcttgttaaca",
  "gtggtcattcaa"
]

# helper functions
def all_equal(arr):
  return arr.count(arr[0]) == len(arr)

def read_fasta(fasta_path):
  record = SeqIO.read(fasta_path, "fasta")
  return record.seq.lower()

def getKmers(sequence, k, step):    
  for x in range(0, len(sequence) - k, step):
    yield sequence[x:x+k]

def makeWhitelist():
  with open(whitelist_file_path, "r") as whitelist_file:
    for record in SeqIO.parse(whitelist_file, "fasta"):
      for kmer in getKmers(record.seq.lower(), K, 1):
        r.sadd(whitelist_key, str(kmer))

def find_conserved_target_kmers():
  alignment = AlignIO.read(ALIGNMENT_PATH, "clustal")
  sequence_ids = [seq.id for seq in alignment]
  print("alignment of", sequence_ids)
  index_of_target = sequence_ids.index(TARGET_ID)
  alignment_length = alignment.get_alignment_length()
  conserved = [1 if all_equal([seq[i] for seq in alignment]) else 0 for i in range(alignment_length)] 
  for start in range(alignment_length - K):
    if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
      continue 
    kmer = str(alignment[index_of_target][start:start+K].seq).lower()
    if "-" in kmer:
      continue
    n_conserved = sum(conserved[start:start+K])
    print(f"{start}:{start+K}", kmer, n_conserved)
    r.zadd("conserved_target_kmers", {kmer: n_conserved})
  most_conserved_kmer = r.zrevrangebyscore("conserved_target_kmers", 9001, 0, withscores=True, start=0, num=1)[0]
  print(f"the most conserved {K}mer is {most_conserved_kmer[0].decode()} with {int(most_conserved_kmer[1])} bases conserved between {sequence_ids}")

def design_guides():
  ncov_consensus = read_fasta(FASTA_PATH)
  kmers = getKmers(str(ncov_consensus), 30, 1)
  guides, count, nconserved = set(), 0, 0
  for kmer in kmers:
    count = count + 1
    if r.sismember(whitelist_key, kmer):
      continue
    if any([seq in kmer for seq in CONSERVED]):
      nconserved = nconserved + 1
      print(count, kmer, "conserved!")
      guides.add(str(Seq(kmer).reverse_complement().transcribe()))
  print("guides", guides, "percent conserved:", nconserved / count * 100)
  with open(OUTFILE_PATH, "w+") as outfile:
    outfile.write("\n".join(list(guides)))

if __name__ == "__main__":
  # makeWhitelist()
  # design_guides()
  find_conserved_target_kmers()



# non-watson-crick base pair map (cDNA right now) NOTE: this was a combinatorial shit-show and took too long
# WOBBLES = {
#   "a": ["t"],
#   "c": ["g"],
#   "g": ["t", "c"],
#   "t": ["g", "a"],
#   "y": ["g", "a"], # y is c or t, c binds g and t binds g or a (in RNA form)
#   "n": ["a", "c", "t", "g"]
# }
# def wobble(sequence, label1, label2):
#   for seq in itertools.product(*[WOBBLES[char] for char in sequence]):
#     string = "".join(seq)
#     print(label1, ":", sequence, "-->", label2, ":", string)
#     yield string

# def design_guides_with_wobble():
#   # read the coronavirus consensus
#   ncov_consensus = read_fasta(FASTA_PATH)
#   # generate kmers
#   kmers = getKmers(str(ncov_consensus), 30, 1)
#   # create guides for each kmer ... and targets for each guide
#   # if a guide has any targets in the whitelist, we skip it
#   i = 0
#   guides = set()
#   for kmer in kmers:
#     print(i, kmer)
#     i = i + 1
#     if r.sismember(whitelist_key, kmer):
#       continue
#     for guide in wobble(kmer, "kmer", "guide"):
#       try:
#         targets = set() 
#         for target in wobble(guide, "guide", "target"):
#           if r.sismember(whitelist_key, target):
#             print("HIT!")
#             raise Exception("guide target in whitelist")
#           targets.add(target)
#         guides.add(guide)
#         print("guide:", guide, "targets:", targets)
#       except Exception as e:
#         print(e)
#         pass
#   print("guides:", guides)
#   with open(os.path.join("guides", target_fasta_file_name), "w") as outfile:
#     outfile.write("\n".join(list(guides)))  

# an even slower way to make wobbles
# def find_wobbles(sequence, nucleotides=['g', 't', 'y', 'n']):
#   return [(i, sequence[i]) for i in range(len(sequence)) if sequence[i] in nucleotides]

# save the whitelist for later
# with open(os.path.join("whitelist", f"lung-cds-{K}mers.txt"), "w") as outfile:
#   outfile.write("\n".join(r.smembers('whitelist:lung:cds')))

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