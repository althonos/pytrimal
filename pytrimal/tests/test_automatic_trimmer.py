import os
import pathlib
import unittest

from .. import Alignment, AutomaticTrimmer



class TestAutomaticTrimmer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data_folder = pathlib.Path(__file__).parent.joinpath("data")

    def _test_method(self, name):
        ali = Alignment.load(self.data_folder.joinpath("ENOG411BWBU.fasta"))
        expected = Alignment.load(self.data_folder.joinpath(f"ENOG411BWBU.{name}.fasta"))

        trimmer = AutomaticTrimmer(name)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_strict_method(self):
        self._test_method("strict")

    def test_strictplus_method(self):
        self._test_method("strictplus")

    def test_automatic1_method(self):
        self._test_method("automated1")

    def test_noallgaps_method(self):
        self._test_method("noallgaps")

    def test_gappyout_method(self):
        self._test_method("gappyout")
