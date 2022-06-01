import os
import pathlib
import unittest

from .. import Alignment, ManualTrimmer


class TestManualTrimmer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data_folder = pathlib.Path(__file__).parent.joinpath("data")

    def _test_parameters(self, gt, cons):
        ali = Alignment.load(self.data_folder.joinpath("ENOG411BWBU.fasta"))
        expected = Alignment.load(self.data_folder.joinpath("ENOG411BWBU.cons{:02}.gt{:02}.fasta".format(cons, int(gt*100))))

        trimmer = ManualTrimmer(gap_threshold=gt, conservation_percentage=cons)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_gap_threshold(self):
        self._test_parameters(gt=0.9, cons=60)
        self._test_parameters(gt=0.4, cons=40)
