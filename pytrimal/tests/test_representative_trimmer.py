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

from .. import _trimal, Alignment, RepresentativeTrimmer


class TestRepresentativeTrimmer(unittest.TestCase):
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

    def _test_representative(self, clusters=None, identity_threshold=None):
        ali = self._load_alignment("ENOG411BWBU.fasta")

        if clusters is not None:
            expected = self._load_alignment(
                "ENOG411BWBU.clusters{}.fasta".format(clusters)
            )
        else:
            expected = self._load_alignment(
                "ENOG411BWBU.maxidentity{}.fasta".format(int(identity_threshold * 100))
            )

        trimmer = RepresentativeTrimmer(
            clusters=clusters, identity_threshold=identity_threshold, backend=self.backend
        )
        trimmed = trimmer.trim(ali)

        if clusters is not None:
            self.assertEqual(len(trimmed.sequences), clusters)
        self.assertTrimmedAlignmentEqual(trimmed, expected)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_clusters5(self):
        self._test_representative(clusters=5)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_clusters10(self):
        self._test_representative(clusters=10)

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_identity75(self):
        self._test_representative(identity_threshold=0.75)

    def test_repr(self):
        trimmer = RepresentativeTrimmer(identity_threshold=0.25)
        self.assertEqual(
            repr(trimmer), "RepresentativeTrimmer(identity_threshold=0.25)"
        )
        trimmer = RepresentativeTrimmer(clusters=2)
        self.assertEqual(repr(trimmer), "RepresentativeTrimmer(clusters=2)")
        trimmer = RepresentativeTrimmer(clusters=3, backend=None)
        self.assertEqual(
            repr(trimmer), "RepresentativeTrimmer(clusters=3, backend=None)"
        )

    def test_pickle(self):
        trimmer = RepresentativeTrimmer(clusters=3, backend=self.backend)
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


class TestRepresentativeTrimmerGeneric(TestRepresentativeTrimmer):
    backend = "generic"


@unittest.skipUnless(_trimal._MMX_RUNTIME_SUPPORT, "MMX not available")
class TestRepresentativeTrimmerMMX(TestRepresentativeTrimmer):
    backend = "mmx"


@unittest.skipUnless(_trimal._SSE2_RUNTIME_SUPPORT, "SSE2 not available")
class TestRepresentativeTrimmerSSE(TestRepresentativeTrimmer):
    backend = "sse"


@unittest.skipUnless(_trimal._AVX2_RUNTIME_SUPPORT, "AVX2 not available")
class TestRepresentativeTrimmerAVX(TestRepresentativeTrimmer):
    backend = "avx"


@unittest.skipUnless(_trimal._NEON_RUNTIME_SUPPORT, "NEON not available")
class TestRepresentativeTrimmerNEON(TestRepresentativeTrimmer):
    backend = "neon"