import math
import os
import unittest

from .. import SimilarityMatrix


class TestSimilarityMatrix(unittest.TestCase):

    def test_init(self):
        self.assertRaises(NotImplementedError, SimilarityMatrix)

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
