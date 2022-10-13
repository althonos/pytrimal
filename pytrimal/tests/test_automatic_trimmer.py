import os
import pickle
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

from .. import _trimal, Alignment, AutomaticTrimmer, SimilarityMatrix


class TestAutomaticTrimmer(unittest.TestCase):

    backend = None

    def assertTrimmedAlignmentEqual(self, trimmed, expected):
        self.assertEqual(len(trimmed.names), len(expected.names))
        self.assertEqual(len(trimmed.sequences), len(expected.sequences))
        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    @staticmethod
    def _load_alignment(name):
        with importlib_resources.path("pytrimal.tests.data", name) as path:
            return Alignment.load(path)

    def _test_method(self, name):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        expected = self._load_alignment("ENOG411BWBU.{}.fasta".format(name))
        trimmer = AutomaticTrimmer(method=name, backend=self.backend)
        trimmed = trimmer.trim(ali)
        self.assertTrimmedAlignmentEqual(trimmed, expected)

    def test_invalid_method(self):
        self.assertRaises(ValueError, AutomaticTrimmer, method="nonsense")
        self.assertRaises(TypeError, AutomaticTrimmer, method=1)

    def test_repr(self):
        trimmer = AutomaticTrimmer("strict")
        self.assertEqual(repr(trimmer), "AutomaticTrimmer('strict')")
        trimmer = AutomaticTrimmer("automated1")
        self.assertEqual(repr(trimmer), "AutomaticTrimmer('automated1')")
        trimmer = AutomaticTrimmer("noduplicateseqs", backend=None)
        self.assertEqual(
            repr(trimmer), "AutomaticTrimmer('noduplicateseqs', backend=None)"
        )

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
    def test_noduplicateseqs_method(self):
        self._test_method("noduplicateseqs")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_custom_similarity_matrix(self):
        alignment = self._load_alignment("ENOG411BWBU.fasta")
        with importlib_resources.open_binary("pytrimal.tests.data", "pam70.json") as f:
            pam70 = SimilarityMatrix(**json.load(f))

        trimmer = AutomaticTrimmer("strict", backend=self.backend)
        trimmed = trimmer.trim(alignment, pam70)

    def test_invalid_characters(self):
        alignment = Alignment([b"seq1", b"seq2"], ["MKKBO", "MKKAY"])
        trimmer = AutomaticTrimmer(method="strict", backend=self.backend)
        self.assertRaises(ValueError, trimmer.trim, alignment)

    def test_pickle(self):
        trimmer = AutomaticTrimmer(method="automated1", backend=self.backend)
        pickled = pickle.loads(pickle.dumps(trimmer))
        ali = Alignment(
            names=[b"Sp8", b"Sp17", b"Sp10", b"Sp26"],
            sequences=[
                "LG-----------TKSD---NNNNNNNNNNNNNNNNWV----------",
                "APDLLL-IGFLLKTV-ATFG-----------------DTWFQLWQGLD",
                "DPAVL--FVIMLGTI-TKFS-----------------SEWFFAWLGLE",
                "AAALLTYLGLFLGTDYENFA-----------------AAAANAWLGLE",
            ],
        )
        t1 = trimmer.trim(ali)
        t2 = pickled.trim(ali)
        self.assertTrimmedAlignmentEqual(t2, t1)


@unittest.skipUnless(_trimal._SSE2_RUNTIME_SUPPORT, "SSE2 not available")
class TestAutomaticTrimmerSSE(TestAutomaticTrimmer):
    backend = "sse"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestAutomaticTrimmerAVX(TestAutomaticTrimmer):
    backend = "avx"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestAutomaticTrimmerNEON(TestAutomaticTrimmer):
    backend = "neon"