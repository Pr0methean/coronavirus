SHELL=/bin/bash

# conda:
# 	conda info -e | grep -q bio || conda create --name bio python=3.7 && \
# 	conda activate bio && pip install -r requirements.txt

install-redis:
	apt install redis-server

redis:
	redis-server

scylla:
	sudo docker run --net=host --name scylla -d scylladb/scylla --experimental 1 --smp 2 --overprovisioned
	
schema:
	docker exec -it scylla cqlsh -e "create keyspace rna with replication = {'class':'SimpleStrategy', 'replication_factor': 1};" && \
	docker exec -it scylla cqlsh -e "create table rna.trie (pre text, next set<text>, primary key (pre));"

test: redis scylla schema
	pytest

clean:
	-docker kill scylla
	-docker container prune -f
	-redis-server stop
