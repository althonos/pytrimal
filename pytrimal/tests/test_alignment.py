import sys
import os
import unittest

try:
    try:
        import importlib.resources as importlib_resources
    except ImportError:
        import importlib_resources
except ImportError:
    importlib_resources = None

from .. import Alignment, TrimmedAlignment


class TestAlignment(unittest.TestCase):

    def setUp(self):
        self.alignment = Alignment(
            names=[b"Sp8", b"Sp10", b"Sp26", b"Sp6", b"Sp17", b"Sp33"],
            sequences=[
                "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
                "-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII",
                "AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI",
                "-----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII",
                "--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF",
                "-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL",
            ],
        )

    def test_load_errors(self):
        self.assertRaises(FileNotFoundError, Alignment.load, "nothing")
        self.assertRaises(IsADirectoryError, Alignment.load, os.getcwd())

    def test_residues(self):
        self.assertEqual(len(self.alignment.residues), 46)
        self.assertEqual(self.alignment.residues[0], "--A---")
        self.assertEqual(self.alignment.residues[10], "IVLLLL")
        with self.assertRaises(IndexError):
            self.alignment.residues[100]
        with self.assertRaises(IndexError):
            self.alignment.residues[46]
        with self.assertRaises(IndexError):
            self.alignment.residues[-100]

    def test_sequences(self):
        self.assertEqual(len(self.alignment.sequences), 6)
        self.assertEqual(self.alignment.sequences[0], "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII")
        self.assertEqual(self.alignment.sequences[4], "--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF")
        self.assertEqual(self.alignment.sequences[-1], "-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL")
        with self.assertRaises(IndexError):
            self.alignment.sequences[100]
        with self.assertRaises(IndexError):
            self.alignment.sequences[6]
        with self.assertRaises(IndexError):
            self.alignment.sequences[-100]

class TestTrimmedAlignment(unittest.TestCase):

    def setUp(self):
        residues_mask = [True] * 46
        residues_mask[:5] = [False]*5
        residues_mask[26:28] = [False]*2
        sequences_mask = [True, True, False, True, True, True]
        self.trimmed = TrimmedAlignment(
            names=[b"Sp8", b"Sp10", b"Sp26", b"Sp6", b"Sp17", b"Sp33"],
            sequences=[
                "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
                "-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII",
                "AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI",
                "-----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII",
                "--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF",
                "-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL",
            ],
            sequences_mask=sequences_mask,
            residues_mask=residues_mask,
        )

    @unittest.skipIf(sys.version_info < (3, 6), "No pathlib support in Python 3.5")
    @unittest.skipUnless(importlib_resources, "importlib.resources not available")
    def test_load(self):
        with importlib_resources.path("pytrimal.tests.data", "example.001.AA.clw") as path:
            trimmed = TrimmedAlignment.load(path)
        self.assertEqual(len(trimmed.sequences), 6)
        self.assertEqual(len(trimmed.residues), 46)
        self.assertEqual(len(trimmed.sequences_mask), 6)
        self.assertEqual(len(trimmed.residues_mask), 46)
        self.assertTrue(all(trimmed.sequences_mask))
        self.assertTrue(all(trimmed.residues_mask))

    def test_original_alignment(self):
        original = self.trimmed.original_alignment()
        self.assertEqual(len(original.sequences), 6)
        self.assertEqual(len(original.residues), 46)

    def test_residues(self):
        self.assertEqual(len(self.trimmed.residues), 39)
        self.assertEqual(self.trimmed.residues[0], "G-AT-")
        with self.assertRaises(IndexError):
            self.trimmed.residues[100]
        with self.assertRaises(IndexError):
            self.trimmed.residues[39]
        with self.assertRaises(IndexError):
            self.trimmed.residues[-100]

    def test_sequences(self):
        self.assertEqual(len(self.trimmed.sequences), 5)
        self.assertEqual(self.trimmed.sequences[3], "TAPDLL-LIGFLLKTVA-TFGDTWFQLWQGLDLNKMPVF")
        self.assertEqual(self.trimmed.sequences[-1], "--PTILNIAGLHMETDI-NFSLAWFQAWGGLEINKQAIL")
        with self.assertRaises(IndexError):
            self.trimmed.sequences[5]
        with self.assertRaises(IndexError):
            self.trimmed.sequences[100]
        with self.assertRaises(IndexError):
            self.trimmed.sequences[-100]
