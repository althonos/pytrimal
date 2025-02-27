import os
import pickle
import sys
import json
import unittest

from .. import _trimal, Alignment, AutomaticTrimmer, SimilarityMatrix
from ._base import TrimmerTestCase, files


class TestAutomaticTrimmer(TrimmerTestCase, unittest.TestCase):
    platform = None

    def _test_method(self, name):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        expected = self._load_alignment("ENOG411BWBU.{}.fasta".format(name))
        trimmer = AutomaticTrimmer(method=name, platform=self.platform)
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
        trimmer = AutomaticTrimmer("noduplicateseqs", platform=None)
        self.assertEqual(
            repr(trimmer), "AutomaticTrimmer('noduplicateseqs', platform=None)"
        )

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_strict_method(self):
        self._test_method("strict")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_strictplus_method(self):
        self._test_method("strictplus")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_automatic1_method(self):
        self._test_method("automated1")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_noallgaps_method(self):
        self._test_method("noallgaps")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_gappyout_method(self):
        self._test_method("gappyout")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_noduplicateseqs_method(self):
        self._test_method("noduplicateseqs")

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_custom_similarity_matrix(self):
        alignment = self._load_alignment("ENOG411BWBU.fasta")
        with self._open_data("pam70.json", "rb") as f:
            pam70 = SimilarityMatrix(**json.load(f))

        trimmer = AutomaticTrimmer("strict", platform=self.platform)
        trimmed = trimmer.trim(alignment, pam70)

    def test_invalid_characters(self):
        alignment = Alignment([b"seq1", b"seq2"], ["MKKBO", "MKKAY"])
        trimmer = AutomaticTrimmer(method="strict", platform=self.platform)
        self.assertRaises(ValueError, trimmer.trim, alignment)

    def test_pickle(self):
        trimmer = AutomaticTrimmer(method="automated1", platform=self.platform)
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
class TestAutomaticTrimmerSSE2(TestAutomaticTrimmer):
    platform = "sse2"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestAutomaticTrimmerAVX2(TestAutomaticTrimmer):
    platform = "avx2"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestAutomaticTrimmerNEON(TestAutomaticTrimmer):
    platform = "neon"