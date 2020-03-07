# Hack Coronavirus with Redis and CRISPR

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Problem: Silent spread of a disease

There is active silent spread of SARS-Cov-2 causing Covid-19 which is 200-400x more lethal than flu.
Silent spread of a disease increases the probability of further increases in virulence / lethality because more hosts enable more replication events and each replication event can yield mutations which cause increased virulence / lethality

## Challenge

How do you design CRISPR guides to target SARS-Cov-2 conserved regions (which are unlikely to mutate) without off-target effects upon the human genome or common probiotic / beneficial bugs.

## Opportunity

design anti-coronavirus CRISPR-Cas13 guides

## Prepare

install clustal omega 
install redis-cli and docker 
download human transcriptome with this link: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz
found on page: https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens
download sequences of SARS, MERS, HKU1, and as many SARS-nCoV-2 genomes as possible from ncbi
generate a consensus between as many SARS-nCoV-2 sequences as can be found at https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/
align the prior outbreak sequences with the new consensus with `clustal -i combined.fasta -o combined.clu -outfmt=clu`
open the alignment and use CTRL-F and type `************` until you can't find any matches, then back up and copy this conserved substring from the SARS-nCoV-2 sequence into some file (i did this already) ... gradually delete stars to get more conserved substrings
run redis for the whitelist `docker run -p 6379:6379 -d redis redis-server --appendonly yes`
combine the files with `cat $(ls -t) > combined.fasta`

## Act

run `python design_guides.py` with python (I'm using 3.8)

## Reflect

fancy stuff didn't work (adding the whole transcriptome to the whitelist was WAY too slow on my laptop...but it's embarassingly parallelizable problem to divide this text file and make kmers)
likewise I tried a wobble powerset algorithm (see comments at bottom of design_guides.py) but this caused a combinatorial explosion because there are an insane number of possible 30mers. Could use a MAX ENTROPY heuristic for this algo to just choose the nucleotides which accomodate most wobble
there are probably some gotchas in the gRNA design which 
since I couldn't whitelist the whole transcriptome, the lung microRNA and long non-coding RNA aren't in the whitelist, which presents an opportunity for off target effects
likewise, I'm not using Levenshtein distance for the whitelist cache hits, only exact matches, and there weren't any. That might be a false negative because guides could tolerate some mispairing, and wobble pairing, but I wanted to produce results and then improve it later

STATUS: Work in Progress Minimum Viable Product

## References

BLASTN Suite

Genbank

Clustal Omega

Tips and Tricks for Cas13 https://zlab.bio/cas13

Freije CA, et al. Programmable inhibition and detection of RNA viruses using Cas13. Mol Cell. 2019. https://doi.org/10.1016/j.molcel.2019.09.013.

pC0046-EF1a-PspCas13b-NES-HIV https://www.addgene.org/103862/

Stemmer, M., Thumberger, T., del Sol Keyer, M., Wittbrodt, J. and Mateo, J.L. CCTop: an intuitive, flexible and reliable CRISPR/Cas9 target prediction tool. PLOS ONE (2015). doi: 10.1371/journal.pone.0124633

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.

Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
