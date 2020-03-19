from unittest import TestCase

from design_guides import conserved_in_alignment, count_conserved
from Bio import AlignIO, SeqIO
from datetime import datetime
from Bio.Seq import Seq

class Test(TestCase):
    alignment = [SeqIO.SeqRecord(i) for i in [
        'ATTAAAGGTTTATCCCTTCCCAGGTAGCAAACCACCCAACTGTCGATCTCTTGTAGGTCTGTCCTCTAAA',
        'CGAACTTGAAAATCTGTGTGCAGGTAGCTCGGCTCCATGCTGAGTGCACTCACGCAGTATAACTAATAAC',
        'TAATTACGGTCGTCGACAGGCAGGTAGTAACTCGCCTATCTGCTGCAGGCTGCTTAGGGTTTCGTCCGTG',
        'TTGCAGCGGATCACCAGCACCAGGTAGTTTCGTCCGGGTGTGACCGAAAGGTAAGAGGGAGACCCTTGTC',
    ]]
    conserved = [int(i) for i in list(
        '0000000100000100000011111110000000100000110000000000000010000010000000')]

    def test_make_hosts(self):
        self.fail()

    def test_make_targets(self):
        self.fail()

    def test_predict_side_effects(self):
        self.fail()

    def test_make_plasmids(self):
        self.fail()

    def test_conserved_in_alignment(self):
        self.assertEqual(conserved_in_alignment(self.alignment, len(self.alignment[0])), self.conserved)

    def test_count_conserved(self):
        # Does not have bases 15-21 conserved
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, 0), ("", 0))
        # Has bases 15-21 conserved
        self.assertEqual(count_conserved(self.alignment, self.conserved, 0, 6),
                         ("ggtttatcccttcccaggtagcaaacca", 9))