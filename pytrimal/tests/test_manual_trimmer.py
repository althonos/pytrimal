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

    def _test_parameters(self, gt, cons):
        with importlib_resources.path("pytrimal.tests.data", "ENOG411BWBU.fasta") as path:
            ali = Alignment.load(path)
        filename = "ENOG411BWBU.cons{:02}.gt{:02}.fasta".format(cons, int(gt*100))
        with importlib_resources.path("pytrimal.tests.data", filename) as path:
            expected = Alignment.load(path)

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
