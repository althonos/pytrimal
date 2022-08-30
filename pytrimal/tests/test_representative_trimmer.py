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

from .. import Alignment, RepresentativeTrimmer
from .._trimal import _SSE2_RUNTIME_SUPPORT


class TestRepresentativeTrimmer(unittest.TestCase):

    @staticmethod
    def _load_alignment(name):
        with importlib_resources.path("pytrimal.tests.data", name) as path:
            return Alignment.load(path)

    def _test_representative(self, clusters=None, identity_threshold=None):
        ali = self._load_alignment("ENOG411BWBU.fasta")

        if clusters is not None:
            expected = self._load_alignment("ENOG411BWBU.clusters{}.fasta".format(clusters))
        else:
            expected = self._load_alignment("ENOG411BWBU.maxidentity{}.fasta".format(int(identity_threshold*100)))

        trimmer = RepresentativeTrimmer(clusters=clusters, identity_threshold=identity_threshold)
        trimmed = trimmer.trim(ali)

        if clusters is not None:
            self.assertEqual(len(trimmed.sequences), clusters)
        self.assertEqual(len(trimmed.names), len(expected.names))
        self.assertEqual(len(trimmed.sequences), len(expected.sequences))
        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

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
        self.assertEqual(repr(trimmer), "RepresentativeTrimmer(identity_threshold=0.25)")
        trimmer = RepresentativeTrimmer(clusters=2)
        self.assertEqual(repr(trimmer), "RepresentativeTrimmer(clusters=2)")
        trimmer = RepresentativeTrimmer(clusters=3, backend=None)
        self.assertEqual(repr(trimmer), "RepresentativeTrimmer(clusters=3, backend=None)")
