import os
import sys
import json
import pickle
import unittest

try:
    try:
        import importlib.resources as importlib_resources
    except ImportError:
        import importlib_resources
except ImportError:
    importlib_resources = None

from .. import _trimal, Alignment, OverlapTrimmer


class TestOverlapTrimmer(unittest.TestCase):
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

    def _test_overlap(self, seq, res):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        expected = self._load_alignment(
            "ENOG411BWBU.seq{}.res{}.fasta".format(seq, res)
        )
        trimmer = OverlapTrimmer(sequence_overlap=seq, residue_overlap=res / 100, backend=self.backend)
        trimmed = trimmer.trim(ali)
        self.assertTrimmedAlignmentEqual(trimmed, expected)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_seqoverlap80_resoverlap80(self):
        self._test_overlap(80, 80)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_seqoverlap40_resoverlap60(self):
        self._test_overlap(40, 60)

    def test_repr(self):
        trimmer = OverlapTrimmer(80, 0.5)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(80.0, 0.5)")
        trimmer = OverlapTrimmer(50, 1.0)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(50.0, 1.0)")
        trimmer = OverlapTrimmer(30, 0.25, backend=None)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(30.0, 0.25, backend=None)")

    def test_pickle(self):
        trimmer = OverlapTrimmer(40, 0.5, backend=self.backend)
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


class TestOverlapTrimmerGeneric(TestOverlapTrimmer):
    backend = "generic"


@unittest.skipUnless(_trimal._SSE2_RUNTIME_SUPPORT, "SSE2 not available")
class TestOverlapTrimmerSSE(TestOverlapTrimmer):
    backend = "sse"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestOverlapTrimmerAVX(TestOverlapTrimmer):
    backend = "avx"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestOverlapTrimmerNEON(TestOverlapTrimmer):
    backend = "neon"