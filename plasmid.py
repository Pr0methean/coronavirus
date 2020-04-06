# why?: therapy requires delivery
# how?: AAV because it's well-researched and established
# what?: pCMV-CasRx + pU6-gRNA

from Bio import SeqIO
import os

# path to plasmid backbone genbank
# BACKBONE_PATH = os.path.join("data", "parts", "AAV-pCMV-EGFP.gb")
# start of EGFP (one-indexed!)
# EGFP_START, EGFP_END = 788, 1507
# path to 5 prime ITR
FIVE_PRIME_ITR_PATH = os.path.join("data", "parts", "AAV-5prime-itr.fa")
# path to CMV promoter/enhancer
PROMOTER_PATH = os.path.join("data", "parts", "pCMV.fa")
# path to kozak ribosomal binding site
KOZAK_PATH = os.path.join("data", "parts", "kozak.fa")
# path to CRISPR DNA
CRISPR_PATH = os.path.join("data", "parts", "crispr", "casrx-dna.fa")
# path to polyA tail
POLYA_PATH = os.path.join("data", "parts", "sv40-poly-a.fa")
# path to 3' ITR
THREE_PRIME_ITR_PATH = os.path.join("data", "parts", "AAV-3prime-itr.fa")
# list of paths
PATHS = [
    FIVE_PRIME_ITR_PATH,
    PROMOTER_PATH,
    KOZAK_PATH,
    CRISPR_PATH,
    POLYA_PATH,
    THREE_PRIME_ITR_PATH
]

def splice(prior, insert, start, end):
    s, e = start - 1, end - 1
    prior_left, removed, prior_right = prior[:s], prior[s:e], prior[e:]
    assert len(prior_left) + len(removed) + len(prior_right) == len(prior)
    return prior_left + insert + prior_right

def add_crispr():
    backbone = SeqIO.read(BACKBONE_PATH, "genbank").seq.lower()
    print(backbone)
    crispr = SeqIO.read(CRISPR_PATH, "fasta").seq.lower()
    print(crispr)
    spliced = splice(backbone, crispr, EGFP_START, EGFP_END)
    print(spliced.index(crispr))
    return spliced

def load_fasta(path, fmt="fasta"):
    return SeqIO.read(path, fmt).seq.lower()

if __name__ == "__main__":
    sequences = [load_fasta(p) for p in PATHS]
    insert_size = sum(map(len, sequences))
    print(insert_size)

