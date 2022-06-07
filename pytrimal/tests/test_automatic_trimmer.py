import os
import sys
import json
import unittest

try:
    try:
        import importlib.resources as importlib_resources
    except ImportError:
        import importlib_resources
except ImportError:
    importlib_resources = None

from .. import Alignment, AutomaticTrimmer, SimilarityMatrix
from .._trimal import _SSE2_RUNTIME_SUPPORT


class TestAutomaticTrimmer(unittest.TestCase):

    @staticmethod
    def _load_alignment(name):
        with importlib_resources.path("pytrimal.tests.data", name) as path:
            return Alignment.load(path)

    def _test_method(self, name):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        expected = self._load_alignment("ENOG411BWBU.{}.fasta".format(name))

        trimmer = AutomaticTrimmer(method=name)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_invalid_method(self):
        self.assertRaises(ValueError, AutomaticTrimmer, method="nonsense")
        self.assertRaises(TypeError, AutomaticTrimmer, method=1)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_strict_method(self):
        self._test_method("strict")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_strictplus_method(self):
        self._test_method("strictplus")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_automatic1_method(self):
        self._test_method("automated1")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_noallgaps_method(self):
        self._test_method("noallgaps")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_gappyout_method(self):
        self._test_method("gappyout")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_custom_similarity_matrix(self):
        alignment = self._load_alignment("ENOG411BWBU.fasta")
        with importlib_resources.open_binary("pytrimal.tests.data", "pam70.json") as f:
            pam70 = SimilarityMatrix(**json.load(f))

        trimmer = AutomaticTrimmer("strict")
        trimmed = trimmer.trim(alignment, pam70)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    @unittest.skipUnless(_SSE2_RUNTIME_SUPPORT, "SSE2 not available")
    def test_consistency_sse2(self):
        ali = self._load_alignment("ENOG411BWBU.fasta")

        trimmer_generic = AutomaticTrimmer(method="automated1", backend=None)
        trimmer_sse = AutomaticTrimmer(method="automated1", backend="sse")

        trimmed_generic = trimmer_generic.trim(ali)
        trimmed_sse = trimmer_sse.trim(ali)

        self.assertEqual(trimmed_generic.names, trimmed_sse.names)
        for seq1, seq2 in zip(trimmed_generic.sequences, trimmed_sse.sequences):
            self.assertEqual(seq1, seq2)

    def test_invalid_character_generic(self):
        alignment = Alignment([b"seq1", b"seq2"], ["MKKBO", "MKKAY"])
        trimmer = AutomaticTrimmer(method="strict", backend=None)
        self.assertRaises(ValueError, trimmer.trim, alignment)

    @unittest.skipUnless(_SSE2_RUNTIME_SUPPORT, "SSE2 not available")
    def test_invalid_character_sse2(self):
        alignment = Alignment([b"seq1", b"seq2"], ["MKKBO", "MKKAY"])
        trimmer = AutomaticTrimmer(method="strict", backend="sse")
        self.assertRaises(ValueError, trimmer.trim, alignment)
