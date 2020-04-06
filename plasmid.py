# why?: therapy requires delivery
# how?: AAV because it's well-researched and established
# what?: pCMV-CasRx + pU6-gRNA

from Bio import SeqIO
import os

# path to plasmid backbone genbank
BACKBONE_PATH = os.path.join("data", "parts", "AAV-pCMV-EGFP.gb")
# end of 3' itr (one-indexed!)
THREE_PRIME_ITR_END = 1926
# path to 5 prime ITR
FIVE_PRIME_ITR_PATH = os.path.join("data", "parts", "AAV-5prime-itr.fa")
# path to CMV promoter/enhancer
PCMV_PATH = os.path.join("data", "parts", "pCMV.fa")
# path to kozak ribosomal binding site
KOZAK_PATH = os.path.join("data", "parts", "kozak.fa")
# path to CRISPR DNA
CRISPR_PATH = os.path.join("data", "parts", "crispr", "casrx-dna.fa")
# path to polyA tail
POLYA_PATH = os.path.join("data", "parts", "sv40-poly-a.fa")
# path to pU6 promoter
PU6_PATH = os.path.join("data", "parts", "pU6.fa")
# path to direct repeat (DR)
DR_PATH = os.path.join("data", "parts", "crispr", "casrx-dr.fa")
# path to 3' ITR
THREE_PRIME_ITR_PATH = os.path.join("data", "parts", "AAV-3prime-itr.fa")
# list of paths
PATHS = [
    FIVE_PRIME_ITR_PATH,
    PCMV_PATH,
    KOZAK_PATH,
    CRISPR_PATH,
    POLYA_PATH,
    THREE_PRIME_ITR_PATH,
    PU6_PATH
]
GUIDE_FILE = "trie_guides_26_full_transcriptome.csv"
GUIDES = [l.rstrip() for l in open(os.path.join("data", "guides", GUIDE_FILE))]
K = len(GUIDES[0])

def splice(prior, insert, start, end):
    s, e = start - 1, end - 1
    prior_left, removed, prior_right = prior[:s], prior[s:e], prior[e:]
    assert len(prior_left) + len(removed) + len(prior_right) == len(prior)
    return prior_left + insert + prior_right

def design_plasmid():
    insert = "".join([str(load_fasta(p)) for p in PATHS])
    dr = load_fasta(DR_PATH)
    g = 0
    while len(insert) < 4700 - len(dr) - K:
        insert = insert + dr + GUIDES[g]
        g = g + 1
    right_side = SeqIO.read(BACKBONE_PATH, "genbank").seq.lower()[THREE_PRIME_ITR_END:]
    plasmid_sequence = insert + right_side
    return plasmid_sequence, len(insert), g

def load_fasta(path, fmt="fasta"):
    return SeqIO.read(path, fmt).seq.lower()

if __name__ == "__main__":
    plasmid_sequence, insert_size, g = design_plasmid()
    print(f">AAV-pCMV-CasRx_w_{g}x{K}_anti-SARS-CoV-2_gRNA_insert_size_{insert_size}")
    print(plasmid_sequence)
    

