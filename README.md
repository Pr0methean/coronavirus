# Hack Coronavirus with Redis and CRISPR

## Problem: Silent spread of a disease

There is active silent spread of SARS-Cov-2 causing Covid-19 which is 100-800x more lethal than flu.
Silent spread of a disease increases the probability of further increases in virulence / lethality because more hosts enable more replication events and each replication event can yield mutations which cause increased virulence / lethality

## Challenge

How do you design CRISPR guides to target SARS-Cov-2 conserved regions (which are unlikely to mutate) without off-target effects upon the human genome or common probiotic / beneficial bugs.

## Opportunity

design anti-coronavirus CRISPR-Cas13 guides

## Prepare

download human transcriptome with this link: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz
found on page: https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens
download sequences of SARS, MERS, HKU1, and SARS-nCov-2 from ncbi
combine the files with `cat $(ls -t) > combined.fasta`
install clustal omega and align the viruses with `clustal -i combined.fasta -o combined.clu -outfmt=clu`
install redis-cli and docker and run `docker run --name some-redis -d redis redis-server --appendonly yes`

## Act

## Reflect

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
