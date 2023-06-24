import os
import pickle
import sys
import unittest

from .. import _trimal, Alignment, ManualTrimmer
from ._base import TrimmerTestCase, files


class TestManualTrimmer(TrimmerTestCase, unittest.TestCase):

    def _test_parameters(self, gt, cons):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        filename = "ENOG411BWBU.cons{:02}.gt{:02}.fasta".format(cons, int(gt * 100))
        expected = self._load_alignment(filename)
        trimmer = ManualTrimmer(gap_threshold=gt, conservation_percentage=cons, backend=self.backend)
        trimmed = trimmer.trim(ali)
        self.assertTrimmedAlignmentEqual(trimmed, expected)

    def test_invalid_parameters(self):
        self.assertRaises(ValueError, ManualTrimmer, gap_threshold=100)
        self.assertRaises(ValueError, ManualTrimmer, gap_threshold=-1)
        self.assertRaises(ValueError, ManualTrimmer, gap_absolute_threshold=-1)
        self.assertRaises(ValueError, ManualTrimmer, conservation_percentage=1000)
        self.assertRaises(ValueError, ManualTrimmer, conservation_percentage=-2)
        self.assertRaises(
            ValueError, ManualTrimmer, gap_threshold=0.5, gap_absolute_threshold=0.5
        )
        self.assertRaises(ValueError, ManualTrimmer, window=5, gap_window=5)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_gap_threshold(self):
        self._test_parameters(gt=0.9, cons=60)
        self._test_parameters(gt=0.4, cons=40)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_window(self):
        ali = self._load_alignment("example.001.AA.clw", "clustal")
        expected = self._load_alignment("example.001.gt90.w3.clw", "clustal")

        trimmer = ManualTrimmer(gap_threshold=0.9, window=3, backend=self.backend)
        trimmed = trimmer.trim(ali)

        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    def test_large_window(self):
        ali = Alignment([b"seq1", b"seq2"], ["M-KKV", "MY-KV"])
        trimmer = ManualTrimmer(gap_threshold=0.9, window=100, backend=self.backend)
        self.assertRaises(Exception, trimmer.trim, ali)

    def test_duplicate_window(self):
        self.assertRaises(ValueError, ManualTrimmer, window=3, gap_window=3)
        self.assertRaises(
            ValueError, ManualTrimmer, window=3, gap_window=3, similarity_window=3
        )

    def test_repr(self):
        trimmer = ManualTrimmer(gap_threshold=0.5)
        self.assertEqual(repr(trimmer), "ManualTrimmer(gap_threshold=0.5)")
        trimmer = ManualTrimmer(window=5, backend=None)
        self.assertEqual(repr(trimmer), "ManualTrimmer(window=5, backend=None)")
        trimmer = ManualTrimmer(
            gap_absolute_threshold=10,
            similarity_threshold=0.5,
            conservation_percentage=50.0,
            gap_window=5,
            similarity_window=5,
            backend=None,
        )
        self.assertEqual(
            repr(trimmer),
            "ManualTrimmer(gap_absolute_threshold=10, similarity_threshold=0.5, conservation_percentage=50.0, gap_window=5, similarity_window=5, backend=None)",
        )

    def test_pickle(self):
        trimmer = ManualTrimmer(gap_threshold=0.4, window=5, backend=self.backend)
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


class TestManualTrimmerGeneric(TestManualTrimmer):
    backend = "generic"


@unittest.skipUnless(_trimal._MMX_RUNTIME_SUPPORT, "MMX not available")
class TestManualTrimmerMMX(TestManualTrimmer):
    backend = "mmx"


@unittest.skipUnless(_trimal._SSE2_RUNTIME_SUPPORT, "SSE2 not available")
class TestManualTrimmerSSE(TestManualTrimmer):
    backend = "sse"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestManualTrimmerAVX(TestManualTrimmer):
    backend = "avx"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestManualTrimmerNEON(TestManualTrimmer):
    backend = "neon"