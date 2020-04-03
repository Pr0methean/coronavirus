# why?: reduce immunogenicity of CRISPR + enable open innovation
# how?: computational directed evolution of immuno-safety
# what?: remove epitopes as gently as possible

import csv
import os

from Bio import Seq, SeqIO

# list of amino acids sorted by hydrophobicity
AAs = list(csv.reader((os.path.join("data", "aas.csv"))))
# name of the fasta file with the CRISPR protein sequence
CRISPR_FASTA_NAME = "PspCas13b.fa"


def read_fasta(fasta_path: str) -> Seq:
    record = SeqIO.read(handle=fasta_path, format="fasta")
    return record.seq.lower()


wild_type = read_fasta(
    os.path.join("data", "parts", "crispr", CRISPR_FASTA_NAME))


# predict MHC1 epitopes with MHCFlurry
def mhc1_epitope_predictions(sequence):
    # TODO
    return []


def get_epitopes(sequence):
    return mhc1_epitope_predictions(sequence)


# replace the most immunogenic amino acid
# with the next less hydrophobic one
def score(aa):
    # TODO
    return 0


def mutate(epitope):
    scores = [score(aa) for aa in epitope]
    target_index = scores.index(max(scores))
    target_letter = epitope[target_index]
    new_aa = AAs[AAs.index(target_letter) + 1]["letter"]
    aas = [aa if i is not target_index else new_aa for i,
                                                       aa in enumerate(
        epitope)]
    return "".join(aas)


# replace an epitope with a mutant
def splice(sequence, epitope, mutant):
    start = sequence.index(epitope)
    end = start + len(epitope) + 1
    spliced = sequence[:start] + mutant + sequence[end:]
    assert len(spliced) is len(sequence)
    return spliced


# mutate all epitopes in a sequence
def evolve(sequence, epitopes):
    new_sequence = sequence
    mutations = [mutate(epitope) for epitope in epitopes]
    for epitope, mutation in zip(epitopes, mutations):
        new_sequence = splice(new_sequence, epitope, mutation)
    return new_sequence


if __name__ == "__main__":
    epitopes = get_epitopes(wild_type)
    sequence = wild_type
    g = 1

    while len(epitopes) > 0:
        sequence = evolve(sequence, epitopes)
        epitopes = get_epitopes(sequence)
        print(f"generation {g}: {len(epitopes)} epitopes")
        g = g + 1

    print("final:\n", sequence)
