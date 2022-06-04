import math
import os
import unittest

from .. import SimilarityMatrix


class TestSimilarityMatrix(unittest.TestCase):

    def test_init_nucleotide(self):
        mx = SimilarityMatrix(
            "ATCG",
            [[5, 0, 0, 4],
             [0, 5, 4, 0],
             [0, 4, 5, 0],
             [4, 0, 0, 5]]
        )
        self.assertEqual(mx.similarity('A', 'A'), 5.0)
        self.assertEqual(mx.similarity('A', 'T'), 0.0)
        self.assertEqual(mx.similarity('A', 'G'), 4.0)

    def test_init_wrong_alphabet_size(self):
        self.assertRaises(
            ValueError,
            SimilarityMatrix,
            "ATC",
            [[5, 0, 0, 4],
             [0, 5, 4, 0],
             [0, 4, 5, 0],
             [4, 0, 0, 5]]
        )

    def test_distance_nt(self):
        matrix = SimilarityMatrix.nt()
        self.assertEqual(matrix.distance('A', 'A'), 0.0)
        self.assertGreater(matrix.distance('A', 'T'), 0.0)
        self.assertRaises(ValueError, matrix.distance, '+', ':')
        self.assertRaises(ValueError, matrix.distance, 'nonsense', 'nonsense')

    def test_distance_aa(self):
        matrix = SimilarityMatrix.aa()
        self.assertEqual(matrix.distance('A', 'A'), 0.0)
        self.assertGreater(matrix.distance('A', 'R'), 0.0)
        self.assertRaises(ValueError, matrix.distance, '+', ':')
