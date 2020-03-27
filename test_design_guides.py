from unittest import TestCase

from Bio import SeqIO

from design_guides import conserved_in_alignment, count_conserved, K, index, host_has


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


# noinspection SpellCheckingInspection
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

    def test_make_hosts(self):
        # TODO
        pass

    def test_make_targets(self):
        # TODO
        pass

    def test_predict_side_effects(self):
        # TODO
        pass

    def test_find(self):
        # TODO
        pass

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

    def test_index(self):
        fake_leveldb = FakeLevelDb()
        index('acctg', fake_leveldb)
        self.assertDictEqual(fake_leveldb.my_dict, {
            b'': b'a',
            b'a': b'c',
            b'ac': b'c',
            b'acc': b't',
            b'acct': b'g',
            b'acctg': b'*'})
        for x in range(2):  # Should be idempotent
            index('accgc', fake_leveldb)
            self.assertDictEqual(fake_leveldb.my_dict, {
                b'': b'a',
                b'a': b'c',
                b'ac': b'c',
                b'acc': b'gt',
                b'accg': b'c',
                b'acct': b'g',
                b'accgc': b'*',
                b'acctg': b'*'})

    def test_host_has(self):
        fake_leveldb = FakeLevelDb()
        index('acctg', fake_leveldb)
        index('ggcat', fake_leveldb)
        self.assertTrue(host_has('acctg', ldb=fake_leveldb, max_mismatches=1))
        self.assertTrue(host_has('ccctg', ldb=fake_leveldb, max_mismatches=1))
        self.assertFalse(host_has('ccctt', ldb=fake_leveldb, max_mismatches=1))
        self.assertTrue(host_has('ccctt', ldb=fake_leveldb, max_mismatches=2))
