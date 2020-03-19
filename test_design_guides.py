from unittest import TestCase

from design_guides import conserved_in_alignment


class Test(TestCase):
    def test_make_hosts(self):
        self.fail()

    def test_make_targets(self):
        self.fail()

    def test_predict_side_effects(self):
        self.fail()

    def test_make_plasmids(self):
        self.fail()

    def test_conserved_in_alignment(self):
        alignment = [
                        'ATTAAAGGTTTATCCCTTCCCAGGTAACAAACCACCCAACTGTCGATCTCTTGTAGGTCTGTCCTCTAAA',
                        'CGAACTTGAAAATCTGTGTGCAGGTCACTCGGCTCCATGCTGAGTGCACTCACGCAGTATAACTAATAAC',
                        'TAATTACGGTCGTCGACAGGCAGGTAGTAACTCGCCTATCTGCTGCAGGCTGCTTAGGGTTTCGTCCGTG',
                        'TTGCAGCGGATCACCAGCACCAGGTGGTTTCGTCCGGGTGTGACCGAAAGGTAAGAGGGAGACCCTTGTC',
        ]
        conserved = [int(i) for i in list(
                        '0000000100000100000011111000000000100000110000000000000010000010000000')]
        self.assertEqual(conserved_in_alignment(alignment, len(alignment[0])), conserved)
