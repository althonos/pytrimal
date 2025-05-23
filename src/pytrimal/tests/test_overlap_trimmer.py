import os
import sys
import json
import pickle
import unittest

from .. import _trimal, Alignment, OverlapTrimmer
from ._base import TrimmerTestCase, files


class TestOverlapTrimmer(TrimmerTestCase, unittest.TestCase):
    platform = None

    def _test_overlap(self, seq, res):
        ali = self._load_alignment("ENOG411BWBU.fasta")
        expected = self._load_alignment(
            "ENOG411BWBU.seq{}.res{}.fasta".format(seq, res)
        )
        trimmer = OverlapTrimmer(sequence_overlap=seq, residue_overlap=res / 100, platform=self.platform)
        trimmed = trimmer.trim(ali)
        self.assertTrimmedAlignmentEqual(trimmed, expected)

    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_seqoverlap80_resoverlap80(self):
        self._test_overlap(80, 80)

    @unittest.skipUnless(files, "importlib.resources.files not available")
    def test_seqoverlap40_resoverlap60(self):
        self._test_overlap(40, 60)

    def test_repr(self):
        trimmer = OverlapTrimmer(80, 0.5)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(80.0, 0.5)")
        trimmer = OverlapTrimmer(50, 1.0)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(50.0, 1.0)")
        trimmer = OverlapTrimmer(30, 0.25, platform=None)
        self.assertEqual(repr(trimmer), "OverlapTrimmer(30.0, 0.25, platform=None)")

    def test_pickle(self):
        trimmer = OverlapTrimmer(40, 0.5, platform=self.platform)
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
class TestOverlapTrimmerSSE2(TestOverlapTrimmer):
    platform = "sse2"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestOverlapTrimmerAVX2(TestOverlapTrimmer):
    platform = "avx2"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestOverlapTrimmerNEON(TestOverlapTrimmer):
    platform = "neon"