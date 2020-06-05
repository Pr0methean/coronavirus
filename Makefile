SHELL=/bin/bash

scylla:
	sudo mkdir -p /var/lib/scylla/data /var/lib/scylla/commitlog
	sudo systemctl unmask docker
	sudo systemctl start docker
	docker-compose up -d
	
schema: scylla
	until docker exec -it scylla cqlsh -e "create keyspace rna with replication = {'class':'SimpleStrategy', 'replication_factor': 1};";	do :; done # work around an issue where Scylla isn't ready yet when docker-compose exits
	docker exec -it scylla cqlsh -e "create table rna.trie (pre text, next set<text>, primary key (pre));"
	docker exec -it scylla cqlsh -e "create table rna.hosts (kmer text, primary key (kmer));"
	docker exec -it scylla cqlsh -e "create table rna.targets (target text, n bigint, start bigint, kmer text, score float, host_has boolean, side_effect boolean, primary key ((target), score, n)) with clustering order by (score desc);"

transcriptome:
	cd data/host && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz

test:  
	pytest

clean:
	-docker kill scylla || true
	-docker container prune -f
	sudo rm -rf /var/lib/scylla

