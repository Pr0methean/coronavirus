import itertools
import os
import pickle
import random
from mmap import mmap, ACCESS_READ

import numpy as np
import redis
from Bio import AlignIO, SeqIO
from sklearn.neighbors import BallTree


def bytesu(string):
    return bytes(string, "UTF-8")


# length of the CRISPR guide RNA
K = 28
# path to a fasta file with host sequences to avoid
HOST_TRANSCRIPTOME = "GCF_000001405.39_GRCh38.p13_rna.fna"  # all RNA in human transcriptome
HOST_LUNG_TISSUE = "lung-tissue-gene-cds.fa"  # just lungs
HOST_FILE = HOST_TRANSCRIPTOME
HOST_PATH = os.path.join("host", HOST_FILE)
# ending token for tries
END = bytesu("*")
EMPTY = bytesu('')
# path to pickle / save the index
REBUILD_INDEX = True
USE_EXISTING_VECTORS = False
INDEX_PATH = "trie.pkl"
# path to alignment and id for the sequence to delete
CORONAVIRUS_CONSENSUS = "HKU1+MERS+SARS+nCoV-Consensus.clu"
TARGET_PATH = os.path.join("alignments", CORONAVIRUS_CONSENSUS)
TARGET_ID = "nCoV"
# mismatch cutoff (e.g. how close must two kmers be to bind?)
CUTOFF = 2
# bases 15-21 of Cas13 gRNA don't tolerate mismatch (Wessels et al 2019)
OFFSET_1, OFFSET_2 = 14, 21
# paths to plasmid part fasta files
PROMOTER_PATH = os.path.join("parts", "pdpn_1_promoter.fa")
DR_SEQUENCE_PATH = os.path.join("parts", "crispr", "dr_sequence_bz_short.fa")
TAIL_PATH = os.path.join("parts", "tail.fa")
# path for output
OUTFILE_PATH = os.path.join("guides", "trie_guides.csv")
VECTORS_PATH = os.path.join("tmp", "vectors")
BASE_A = 0
BASE_C = 1
BASE_G = 2
BASE_T = 3
WILDCARD_EXPANSION = {
    ord('a'): {BASE_A},
    ord('c'): {BASE_C},
    ord('g'): {BASE_G},
    ord('t'): {BASE_T},
    ord('w'): {BASE_A, BASE_T},
    ord('n'): {BASE_A, BASE_C, BASE_G, BASE_T}
}
VEC_TO_KMER = {
    BASE_A: 'a',
    BASE_C: 'c',
    BASE_G: 'g',
    BASE_T: 't'
}


def all_equal(arr):
    return arr.count(arr[0]) == len(arr)


def getKmers(sequence: str, k: int, step: int):
    for x in range(0, len(sequence) - k + 1, step):
        yield sequence[x:x + k]


def vec2kmer(vec: bytes):
    return ''.join(VEC_TO_KMER[byte] for byte in vec)


def kmer2vecs(kmer: bytes):
    return (bytes(vec) for vec in itertools.product(*[WILDCARD_EXPANSION[base] for base in kmer]))


def load_from_pickle(file_path):
    with open(file_path, "rb") as pickle_file:
        unpickled = pickle.load(pickle_file)
    return unpickled


def byteses2array(byteses, k):
    num_byteses = len(byteses)
    concatenated = bytearray(num_byteses * k)
    for i in range(0, num_byteses):
        concatenated[i * k:(i + 1) * k] = byteses[i]
    array = np.frombuffer(concatenated, dtype='uint8')
    array = array.reshape(num_byteses, k)
    return array


def host_has(kmer: str, tree: BallTree, max_mismatches=CUTOFF, k=K):
    distance, closest = tree.query(byteses2array(list(kmer2vecs(bytesu(kmer))), k), 1, return_distance=True)
    distance = int(distance * k)
    print(f"closest to {kmer} is {vec2kmer(tree.data[closest[0][0]])} at distance {distance}")
    if distance > max_mismatches:
        print(f"allow")
        return False
    print(f"avoid")
    return True


def make_hosts(input_path=HOST_PATH, k=K, index_path=INDEX_PATH, vectors_path=VECTORS_PATH,
               use_existing_vectors=USE_EXISTING_VECTORS):
    with open(input_path, "r") as host_file:
        if not use_existing_vectors:
            kmers = {bytesu(str(kmer))
                     for record in SeqIO.parse(host_file, "fasta")
                     for kmer in getKmers(record.seq.lower(), k=k, step=1)}
            print(f"{len(kmers)} unique kmers")
            if os.path.exists(vectors_path):
                os.remove(vectors_path)
            with open(vectors_path, "wb+") as temp:
                [temp.write(vec)
                 for kmer in kmers
                 for vec in kmer2vecs(kmer)]
            del kmers
        size = os.path.getsize(vectors_path)
        num_vectors = int(size / k)
        print(f"{num_vectors} unique vectors")
        array = np.frombuffer(open_as_mmap(vectors_path), dtype='uint8')
        print("Array made")
        array = array.reshape(int(num_vectors), int(k))
        tree = BallTree(array, metric='hamming')
        with open(index_path, "wb") as index_file:
            pickle.dump(tree, index_file)
        return tree


def open_as_mmap(file_path):
    temp_file_no = os.open(file_path, os.O_RDONLY)
    mmapped = mmap(temp_file_no, length=os.path.getsize(file_path), access=ACCESS_READ)
    return mmapped


def make_targets(db, target_path=TARGET_PATH, target_id=TARGET_ID, k=K):
    targets_key = f"targets_{k}"
    alignment = AlignIO.read(target_path, "clustal")
    seq_ids = [seq.id for seq in alignment]
    index_of_target = seq_ids.index(target_id)
    alignment_length = alignment.get_alignment_length()
    conserved = conserved_in_alignment(alignment, alignment_length)
    for start in range(alignment_length - k + 1):
        kmer, n_conserved = count_conserved(alignment, conserved, index_of_target, start, k)
        if n_conserved > 0:
            print(f"{kmer} at {start} has {int(n_conserved)} conserved bases")
            db.zadd(targets_key, {kmer: n_conserved})
    most = db.zrevrangebyscore(targets_key, 9001, 0, withscores=True, start=0, num=1)[0]
    print(
        f"the most conserved {K}mer {most[0].decode()} has {int(most[1])} bases conserved in {seq_ids}")


def conserved_in_alignment(alignment, alignment_length):
    return [1 if all_equal(
        [seq[i] for seq in alignment]) else 0 for i in range(alignment_length)]


def count_conserved(alignment, conserved, index_of_target, start, k=K):
    if not all(conserved[start + OFFSET_1:start + OFFSET_2]):
        return "", 0
    else:
        kmer = str(alignment[index_of_target][start:start + k].seq).lower()
        if "-" in kmer:
            return "", 0
        n_conserved = sum(conserved[start:start + k])
    return kmer, n_conserved


def predict_side_effects(tree, db, out_path=OUTFILE_PATH, k=K, max_mismatches=CUTOFF):
    targets_key = f"targets_{k}"
    good_targets_key = f"good_targets_{k}"
    targets = db.zrevrangebyscore(targets_key, 9001, 0)
    for target in targets:
        t = target.decode()
        should_avoid = host_has(t, tree, max_mismatches, k)
        if should_avoid:
            continue
        db.zadd(good_targets_key, {t: db.zscore(targets_key, t)})
    with open(out_path, "w+") as outfile:
        for k, good_target in enumerate(db.zrevrangebyscore(good_targets_key, 90, 0)):
            good_target_string = good_target.decode()
            print("good target", k, good_target_string)
            outfile.write(good_target_string + "\n")
    print(f"saved {db.zcard(good_targets_key)} good targets at {out_path}")


def main(r, rebuild_index=REBUILD_INDEX, host_path=HOST_PATH, target_path=TARGET_PATH,
         target_id=TARGET_ID, k=K, index_path=INDEX_PATH, vectors_path=VECTORS_PATH,
         out_path=OUTFILE_PATH, use_existing_vectors=USE_EXISTING_VECTORS, max_mismatches_for_side_effect=CUTOFF):
    tree = None
    if rebuild_index:
        tree = make_hosts(host_path, k, index_path=index_path, vectors_path=vectors_path,
                          use_existing_vectors=use_existing_vectors)
    else:
        tree = load_from_pickle(index_path)
    # test the trie lookup works
    for i in range(5):
        host_has(vec2kmer(random.choice(tree.data)), tree=tree, k=k,
                 max_mismatches=max_mismatches_for_side_effect)
    make_targets(db=r, target_path=target_path, target_id=target_id, k=k)
    predict_side_effects(db=r, tree=tree, out_path=out_path, k=k,
                         max_mismatches=max_mismatches_for_side_effect)


if __name__ == "__main__":
    main(r=redis.Redis(host='localhost', port=6379))
