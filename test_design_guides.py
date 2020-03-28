import copy
import os
import tempfile
from unittest import TestCase

from Bio import SeqIO

from design_guides import conserved_in_alignment, count_conserved, K, index, host_has, make_hosts, getKmers, find, \
    make_targets, bytesu, predict_side_effects


class FakeWriteBatch:
    def __init__(self, leveldb):
        self.leveldb = leveldb

    def __enter__(self):
        self.my_dict: dict[bytes, bytes] = {}
        return self

    def put(self, key, value):
        self.my_dict[key] = value

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.leveldb.my_dict.update(self.my_dict)
        self.my_dict = None


class FakeLevelDb:
    def __init__(self):
        self.my_dict = {}

    def get(self, key, default):
        if key in self.my_dict:
            return self.my_dict[key]
        else:
            return default

    def write_batch(self, transaction=False) -> FakeWriteBatch:
        wb = FakeWriteBatch(self)
        return wb


class FakeRedis:
    def __init__(self):
        self.my_dict = {}

    def sadd(self, key: str, value: str) -> bool:
        if key not in self.my_dict:
            self.my_dict[key] = {value}
            return True
        else:
            return self.my_dict[key].add(value)

    def exists(self, key: str):
        return 1 if key in self.my_dict else 0

    def zadd(self, key: str, values: dict):
        if key not in self.my_dict:
            self.my_dict[key] = copy.copy(values)
        else:
            self.my_dict[key].update(values)

    def _zrevrangebyscore(self, name: str, max, min, start=None, num=None,
                          withscores=False, score_cast_func=float):
        items_by_score = sorted(self.my_dict[name].items(), key=lambda item: -item[1])
        if start:
            items_by_score = items_by_score[start:]
        if num:
            items_by_score = items_by_score[:num]
        for key, score in items_by_score:
            if min <= score <= max:
                if withscores:
                    yield bytesu(key), score_cast_func(score)
                else:
                    yield bytesu(key)

    def zrevrangebyscore(self, name: str, max, min, start=None, num=None,
                         withscores=False, score_cast_func=float):
        return list(self._zrevrangebyscore(name, max, min, start, num, withscores, score_cast_func))

    def zscore(self, name: str, key: str):
        return self.my_dict[name][key]

    def zcard(self, name: str):
        return len(self.my_dict[name])

    def get(self, name: str):
        if name in self.my_dict:
            return copy.copy(self.my_dict[name])
        else:
            return None


# noinspection SpellCheckingInspection,PyTypeChecker
class Test(TestCase):
    alignment = [SeqIO.SeqRecord(i) for i in [
        'ATTAAAGGTTTATCCCTTCCCAGGTAGCAAACCACCCAACTGTCGATCTCTTGTAGGTCTGTCCTCTAAA',
        'CGAACTTGAAAATCTGTGTGCAGGTAGCTCGGCTCCATGCTGTCGACACTCACGCAGTATAACTAATAAC',
        'TAATTACGGTCGTCGACAGGCAGGTAGTAACTCGCCTATCTGTCGAAGGCTGCTTAGGGTTTCGTCCGTG',
        'TTGCAGCGGATCACCAGCACCAGGTAGTTTCGTCCGGGTGTGTCGAAAAGGTAAGAGGGAGACCCTTGTC',
    ]]
    conserved = [int(i) for i in list(
        '0000000100000100000011111110000000100000111111000000000010000010000000')]
    alignment_length = len(alignment[0])

    def test_get_kmers(self):
        self.assertEqual(list(getKmers('acacaacc', 5, 1)),
                         ['acaca', 'cacaa', 'acaac', 'caacc'])

    def test_conserved_in_alignment(self):
        self.assertEqual(conserved_in_alignment(self.alignment, self.alignment_length), self.conserved)

    def test_count_conserved(self):
        # End cases
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, 0), ("", 0))
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, self.alignment_length - K + 1), ("", 0))
        # Has bases 15-20 but not 21 conserved
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, 26), ("", 0))
        # Has bases 15-21 conserved
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, 6),
                         ("ggtttatcccttcccaggtagcaaacca", 9))
        # Has dashes
        self.assertEqual(
            count_conserved([SeqIO.SeqRecord("ggtttatcccttcccaggtagcaaacc-")], "ggtttatcccttcccaggtagcaaacc-", 0, 0),
            ("", 0))

    def test_index(self):
        fake_redis = FakeRedis()
        index('acctg', fake_redis)
        self.assertDictEqual(fake_redis.my_dict, {
            '': {'a'},
            'a': {'c'},
            'ac': {'c'},
            'acc': {'t'},
            'acct': {'g'},
            'acctg': {'*'}})
        for x in range(2):  # Should be idempotent
            index('accgc', fake_redis)
            self.assertDictEqual({
                '': {'a'},
                'a': {'c'},
                'ac': {'c'},
                'acc': {'g', 't'},
                'accg': {'c'},
                'acct': {'g'},
                'accgc': {'*'},
                'acctg': {'*'},
            }, fake_redis.my_dict)

    def test_find(self):
        fake_redis = FakeRedis()
        index('acctg', fake_redis)
        index('ggcat', fake_redis)
        self.assertEqual(len(list(find('acctg', tdb=fake_redis, max_mismatches=1, k=5))), 1)
        self.assertEqual(len(list(find('acgtg', tdb=fake_redis, max_mismatches=1, k=5))), 1)
        self.assertEqual(len(list(find('acgcg', tdb=fake_redis, max_mismatches=1, k=5))), 0)

    def test_host_has(self):
        fake_redis = FakeRedis()
        index('acctg', fake_redis)
        index('ggcat', fake_redis)
        self.assertTrue(host_has('acctg', tdb=fake_redis, max_mismatches=1, k=5))
        self.assertTrue(host_has('ccctg', tdb=fake_redis, max_mismatches=1, k=5))
        self.assertFalse(host_has('ccctt', tdb=fake_redis, max_mismatches=1, k=5))
        self.assertTrue(host_has('ccctt', tdb=fake_redis, max_mismatches=2, k=5))

    def test_make_hosts(self):
        test_host_path = os.path.join("testdata", "unit_test_host.fa")
        fake_trie_redis = FakeRedis()
        fake_redis = FakeRedis()
        make_hosts(input_path=test_host_path, db=fake_redis, tdb=fake_trie_redis, k=5)
        self.assertEqual({
            '': {'a', 'c', 'g', 't'},
            'a': {'c', 'g'},
            'ac': {'a'},
            'aca': {'a', 'c'},
            'acaa': {'c'},
            'acaac': {'*'},
            'acac': {'a'},
            'acaca': {'*'},
            'ag': {'g'},
            'agg': {'t'},
            'aggt': {'a'},
            'aggta': {'*'},
            'c': {'a'},
            'ca': {'a', 'c'},
            'caa': {'c'},
            'caac': {'c'},
            'caacc': {'*'},
            'cac': {'a'},
            'caca': {'a'},
            'cacaa': {'*'},
            'g': {'g', 't'},
            'gg': {'t'},
            'ggt': {'a'},
            'ggta': {'g'},
            'ggtag': {'*'},
            'gt': {'a'},
            'gta': {'g'},
            'gtag': {'a'},
            'gtaga': {'*'},
            't': {'a'},
            'ta': {'g'},
            'tag': {'g'},
            'tagg': {'t'},
            'taggt': {'*'},
        }, fake_trie_redis.my_dict)
        self.assertEqual({
            'acaca', 'cacaa', 'acaac', 'caacc',
            'taggt', 'aggta', 'ggtag', 'gtaga',
        }, fake_redis.my_dict["hosts"])

    def test_make_targets(self):
        fake_redis = FakeRedis()
        test_target_path = os.path.join("testdata", "unit_test_target.clu")
        make_targets(db=fake_redis, target_path=test_target_path, target_id="nCoV", k=28)
        self.assertEqual([
            (b"caaccaactttcgatctcttggtagatc", 11.0),
            (b"cggtcgtcgacaggcaggtagtaactcg", 9.0),
        ], list(fake_redis.zrevrangebyscore('targets_28', 9001, 0, withscores=True)))

    def test_predict_side_effects(self):
        fake_redis = FakeRedis()
        fake_trie_redis = FakeRedis()
        fake_redis.zadd('targets_5', {
            "caacc": 5.0,
            "cggtc": 4.0
        })
        index("gggtc", fake_trie_redis)
        index("atctg", fake_trie_redis)
        name = None
        with tempfile.NamedTemporaryFile(delete=False) as outfile:
            name = outfile.name
            predict_side_effects(db=fake_redis, out_path=outfile.name, tdb=fake_trie_redis, k=5,
                                 max_mismatches=1)
        with open(name, "r") as outfile:
            self.assertEqual(["caacc\n"], list(outfile.readlines()))
