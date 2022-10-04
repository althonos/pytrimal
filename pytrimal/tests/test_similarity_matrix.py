import math
import os
import sys
import unittest

from .. import SimilarityMatrix


class TestSimilarityMatrix(unittest.TestCase):
    def test_init_nucleotide(self):
        mx = SimilarityMatrix(
            "ATCG", [[5, 0, 0, 4], [0, 5, 4, 0], [0, 4, 5, 0], [4, 0, 0, 5]]
        )
        self.assertEqual(mx.similarity("A", "A"), 5.0)
        self.assertEqual(mx.similarity("A", "T"), 0.0)
        self.assertEqual(mx.similarity("A", "G"), 4.0)

    def test_init_wrong_alphabet_size(self):
        self.assertRaises(
            ValueError,
            SimilarityMatrix,
            "ATC",
            [[5, 0, 0, 4], [0, 5, 4, 0], [0, 4, 5, 0], [4, 0, 0, 5]],
        )

    def test_length(self):
        aa = SimilarityMatrix.aa()
        self.assertEqual(len(aa), 20)
        nt = SimilarityMatrix.nt()
        self.assertEqual(len(nt), 5)
        dn = SimilarityMatrix.nt(degenerated=True)
        self.assertEqual(len(dn), 15)

    def test_distance_nt(self):
        matrix = SimilarityMatrix.nt()
        self.assertEqual(matrix.distance("A", "A"), 0.0)
        self.assertGreater(matrix.distance("A", "T"), 0.0)
        self.assertRaises(ValueError, matrix.distance, "+", ":")
        self.assertRaises(ValueError, matrix.distance, "nonsense", "nonsense")

    def test_distance_aa(self):
        matrix = SimilarityMatrix.aa()
        self.assertEqual(matrix.distance("A", "A"), 0.0)
        self.assertGreater(matrix.distance("A", "R"), 0.0)
        self.assertRaises(ValueError, matrix.distance, "+", ":")

    @unittest.skipUnless(
        sys.implementation.name == "cpython",
        "buffer protocol only supported on CPython",
    )
    def test_memoryview(self):
        matrix = SimilarityMatrix.nt()
        mv = memoryview(matrix)
        self.assertEqual(mv.ndim, 2)
        self.assertEqual(mv.format, "f")
        self.assertTrue(mv.readonly)
        self.assertEqual(mv[0, 0], 1.0)
        self.assertEqual(mv[0, 1], 0.0)
