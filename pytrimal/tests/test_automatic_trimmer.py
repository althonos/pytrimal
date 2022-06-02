import os
import unittest

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

from .. import Alignment, AutomaticTrimmer


class TestAutomaticTrimmer(unittest.TestCase):

    def _test_method(self, name):
        with importlib_resources.path("pytrimal.tests.data", "ENOG411BWBU.fasta") as path:
            ali = Alignment.load(path)
        with importlib_resources.path("pytrimal.tests.data", "ENOG411BWBU.{}.fasta".format(name)) as path:
            expected = Alignment.load(path)

        trimmer = AutomaticTrimmer(name)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_invalid_method(self):
        self.assertRaises(ValueError, AutomaticTrimmer, method="nonsense")
        self.assertRaises(TypeError, AutomaticTrimmer, method=1)

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
