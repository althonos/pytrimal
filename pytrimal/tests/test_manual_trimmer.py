import os
import sys
import unittest

try:
    try:
        import importlib.resources as importlib_resources
    except ImportError:
        import importlib_resources
except ImportError:
    importlib_resources = None

from .. import Alignment, ManualTrimmer


class TestManualTrimmer(unittest.TestCase):

    @staticmethod
    def _load_alignment(name):
        with importlib_resources.path("pytrimal.tests.data", name) as path:
            return Alignment.load(path)

    def _test_parameters(self, gt, cons):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        filename = "ENOG411BWBU.cons{:02}.gt{:02}.fasta".format(cons, int(gt*100))
        expected = self._load_alignment(filename)

        trimmer = ManualTrimmer(gap_threshold=gt, conservation_percentage=cons)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_invalid_parameters(self):
        self.assertRaises(ValueError, ManualTrimmer, gap_threshold=100)
        self.assertRaises(ValueError, ManualTrimmer, gap_threshold=-1)
        self.assertRaises(ValueError, ManualTrimmer, gap_absolute_threshold=-1)
        self.assertRaises(ValueError, ManualTrimmer, conservation_percentage=1000)
        self.assertRaises(ValueError, ManualTrimmer, conservation_percentage=-2)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_gap_threshold(self):
        self._test_parameters(gt=0.9, cons=60)
        self._test_parameters(gt=0.4, cons=40)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_window(self):
        ali = self._load_alignment("example.001.AA.clw")
        expected = self._load_alignment("example.001.gt90.w3.fasta")

        trimmer = ManualTrimmer(gap_threshold=0.9, window=3)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_large_window(self):
        ali = Alignment([b"seq1", b"seq2"], ["M-KKV", "MY-KV"])
        trimmer = ManualTrimmer(gap_threshold=0.9, window=100)
        self.assertRaises(Exception, trimmer.trim, ali)

    def test_duplicate_window(self):
        self.assertRaises(ValueError, ManualTrimmer, window=3, gap_window=3)
        self.assertRaises(ValueError, ManualTrimmer, window=3, gap_window=3, similarity_window=3)
