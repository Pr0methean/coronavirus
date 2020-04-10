SHELL=/bin/bash

install-redis:
	apt install redis-server

redis:
	redis-server

scylla:
	sudo docker run --net=host --name scylla scylladb/scylla --experimental 1 --smp 2 --overprovisioned 1
	
schema:
	docker exec -it scylla cqlsh -e "create keyspace rna with replication = {'class':'SimpleStrategy', 'replication_factor': 1};" && \
	docker exec -it scylla cqlsh -e "create table rna.trie (pre text, next set<text>, primary key (pre));"

data:
	cd data/host && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz

test: redis scylla schema
	pytest

clean:
	-docker kill scylla
	-docker container prune -f
	-redis-server stop
