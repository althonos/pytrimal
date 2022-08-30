import io
import os
import sys
import unittest
import textwrap
import tempfile

try:
    import Bio.Align
    import Bio.SeqRecord
    import Bio.Seq
except ImportError:
    Bio = None

try:
    import pyhmmer.easel
except ImportError:
    pyhmmer = None

from .. import Alignment, TrimmedAlignment


DATA = {
    "clustal": textwrap.dedent(
        """
        CLUSTAL 2.0.12 multiple sequence alignment


        Sp8             -----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII
        Sp10            -------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII
        Sp26            AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI
        Sp6             -----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII
        Sp17            --FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF
        Sp33            -------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL
                                  :    : : *    :.        * **:::    :
        """
    ),
    "fasta": textwrap.dedent(
        """
        >Sp8
        -----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII
        >Sp10
        -------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII
        >Sp26
        AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI
        >Sp6
        -----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII
        >Sp17
        --FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF
        >Sp33
        -------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL
        """
    ),
    "nexus": textwrap.dedent(
        """
        #NEXUS
        BEGIN DATA;
         DIMENSIONS NTAX=6 NCHAR=46;
        FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-;
        [Name: Sp8     Len: 46]
        [Name: Sp10    Len: 46]
        [Name: Sp26    Len: 46]
        [Name: Sp6     Len: 46]
        [Name: Sp17    Len: 46]
        [Name: Sp33    Len: 46]

        MATRIX
        Sp8      -----GLGKV IV-YGIVLGT KSDQFSNWVV WLFPWNGLQI HMMGII
        Sp10     -------DPA VL-FVIMLGT IT-KFS--SE WFFAWLGLEI NMMVII
        Sp26     AAAAAAAAAL LTYLGLFLGT DYENFA--AA AANAWLGLEI NMMAQI
        Sp6      -----ASGAI LT-LGIYLFT LCAVIS--VS WYLAWLGLEI NMMAII
        Sp17     --FAYTAPDL L-LIGFLLKT VA-TFG--DT WFQLWQGLDL NKMPVF
        Sp33     -------PTI LNIAGLHMET DI-NFS--LA WFQAWGGLEI NKQAIL

        ;
        END;

        """
    ),
    "pir": textwrap.dedent(
        """
        >P1;Sp8
        TEST SEQUENCE SP8
          -----GLGKV IV-YGIVLGT KSDQFSNWVV WLFPWNGLQI HMMGII*

        >P1;Sp10
        TEST SEQUENCE SP10
          -------DPA VL-FVIMLGT IT-KFS--SE WFFAWLGLEI NMMVII*

        >P1;Sp26
        TEST SEQUENCE SP26
          AAAAAAAAAL LTYLGLFLGT DYENFA--AA AANAWLGLEI NMMAQI*

        >P1;Sp6
        TEST SEQUENCE SP6
          -----ASGAI LT-LGIYLFT LCAVIS--VS WYLAWLGLEI NMMAII*

        >P1;Sp17
        TEST SEQUENCE SP17
          --FAYTAPDL L-LIGFLLKT VA-TFG--DT WFQLWQGLDL NKMPVF*

        >P1;Sp33
        TEST SEQUENCE SP33
          -------PTI LNIAGLHMET DI-NFS--LA WFQAWGGLEI NKQAIL*

        """
    ),
    "phylip": textwrap.dedent(
        """
         6 46
        Sp8          -----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII
        Sp10         -------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII
        Sp26         AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI
        Sp6          -----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII
        Sp17         --FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF
        Sp33         -------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL


        """
    ),
    "phylip32": textwrap.dedent(
        """
         6 46
        Sp8          -----GLGKV IV-YGIVLGT KSDQFSNWVV WLFPWNGLQI HMMGII

        Sp10         -------DPA VL-FVIMLGT IT-KFS--SE WFFAWLGLEI NMMVII

        Sp26         AAAAAAAAAL LTYLGLFLGT DYENFA--AA AANAWLGLEI NMMAQI

        Sp6          -----ASGAI LT-LGIYLFT LCAVIS--VS WYLAWLGLEI NMMAII

        Sp17         --FAYTAPDL L-LIGFLLKT VA-TFG--DT WFQLWQGLDL NKMPVF

        Sp33         -------PTI LNIAGLHMET DI-NFS--LA WFQAWGGLEI NKQAIL

        """
    ),
}


class TestAlignment(unittest.TestCase):

    type = Alignment

    def setUp(self):
        self.alignment = self.type(
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

    def test_dump_error(self):
        ali = Alignment([b"seq1", b"seq2"], ["MVVK", "MVYK"])
        self.assertRaises(FileNotFoundError, ali.dump, "/some/nonsensical/path")
        self.assertRaises(IsADirectoryError, ali.dump, os.getcwd())
        self.assertRaises(TypeError, ali.dump, io.StringIO())

    def test_dump_fileobj(self):
        ali = Alignment([b"seq1", b"seq2"], ["MVVK", "MVYK"])
        s = io.BytesIO()
        ali.dump(s)
        self.assertEqual(
            s.getvalue().decode().splitlines(), [">seq1", "MVVK", ">seq2", "MVYK"]
        )

    def test_dump_filename(self):
        ali = Alignment([b"seq1", b"seq2"], ["MVVK", "MVYK"])
        s = ali.dumps()
        self.assertEqual(s.splitlines(), [">seq1", "MVVK", ">seq2", "MVYK"])

    def test_dumps(self):
        ali = Alignment([b"seq1", b"seq2"], ["MVVK", "MVYK"])
        s = ali.dumps()
        self.assertEqual(s.splitlines(), [">seq1", "MVVK", ">seq2", "MVYK"])

    def _test_load_filename(self, format):
        with tempfile.NamedTemporaryFile(suffix=format, mode="wb") as tmp:
            tmp.write(DATA[format].lstrip().encode())
            tmp.flush()
            ali = self.type.load(tmp.name)
        self.assertEqual(ali.names, self.alignment.names)
        self.assertEqual(list(ali.sequences), list(self.alignment.sequences))

    def _test_load_fileobj(self, format):
        data = io.BytesIO(DATA[format].lstrip().encode())
        ali = self.type.load(data, format)
        self.assertEqual(ali.names, self.alignment.names)
        self.assertEqual(list(ali.sequences), list(self.alignment.sequences))

    def test_load_filename_fasta(self):
        self._test_load_filename("fasta")

    def test_load_filename_clustal(self):
        self._test_load_filename("clustal")

    def test_load_filename_phylip(self):
        self._test_load_filename("phylip")

    def test_load_filename_phylip32(self):
        self._test_load_filename("phylip32")

    def test_load_filename_nexus(self):
        self._test_load_filename("nexus")

    def test_load_fileobj_fasta(self):
        self._test_load_fileobj("fasta")

    def test_load_fileobj_clustal(self):
        self._test_load_fileobj("clustal")

    def test_load_fileobj_phylip(self):
        self._test_load_fileobj("phylip")

    def test_load_fileobj_pir(self):
        self._test_load_fileobj("pir")

    def test_load_fileobj_nexus(self):
        self._test_load_fileobj("nexus")

    def test_load_errors(self):
        self.assertRaises(FileNotFoundError, self.type.load, "nothing")
        self.assertRaises(IsADirectoryError, self.type.load, os.getcwd())
        self.assertRaises(TypeError, self.type.load, io.StringIO(), "fasta")

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

    def test_residues_slice(self):
        res = self.alignment.residues
        self.assertEqual(list(res[:30:3]), list(res)[:30:3])
        self.assertEqual(list(res[:-1:7]), list(res)[:-1:7])
        self.assertTrue(res[:][:2])

        empty = self.type([], [])
        self.assertFalse(empty.residues[:])
        self.assertFalse(empty.residues[:][:2])

    def test_sequences(self):
        self.assertEqual(len(self.alignment.sequences), 6)
        self.assertEqual(
            self.alignment.sequences[0],
            "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
        )
        self.assertEqual(
            self.alignment.sequences[4],
            "--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF",
        )
        self.assertEqual(
            self.alignment.sequences[-1],
            "-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL",
        )
        with self.assertRaises(IndexError):
            self.alignment.sequences[100]
        with self.assertRaises(IndexError):
            self.alignment.sequences[6]
        with self.assertRaises(IndexError):
            self.alignment.sequences[-100]

    def test_sequences_slice(self):
        seqs = self.alignment.sequences
        self.assertEqual(list(seqs[:5:2]), list(seqs)[:5:2])
        self.assertEqual(list(seqs[:-1:2]), list(seqs)[:-1:2])
        self.assertTrue(seqs[:][:2])

        empty = self.type([], [])
        self.assertFalse(empty.sequences[:])
        self.assertFalse(empty.sequences[:][:2])

    @unittest.skipUnless(Bio, "Biopython is required for this test")
    def test_from_biopython(self):
        SeqRecord = Bio.SeqRecord.SeqRecord
        Seq = Bio.Seq.Seq
        records = [
            SeqRecord(Seq("-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII"), "Sp8"),
            SeqRecord(Seq("-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII"), "Sp10"),
            SeqRecord(Seq("AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI"), "Sp26"),
            SeqRecord(Seq("-----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII"), "Sp6"),
            SeqRecord(Seq("--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF"), "Sp17"),
            SeqRecord(Seq("-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL"), "Sp33"),
        ]
        alignment = Alignment.from_biopython(records)
        self.assertEqual(list(alignment.sequences), list(self.alignment.sequences))
        self.assertEqual(list(alignment.names), list(self.alignment.names))

    @unittest.skipUnless(Bio, "Biopython is required for this test")
    def test_to_biopython(self):
        msa = self.alignment.to_biopython()
        self.assertEqual([str(r.seq) for r in msa], list(self.alignment.sequences))
        self.assertEqual([r.id for r in msa], [name.decode() for name in self.alignment.names])

    @unittest.skipUnless(pyhmmer, "PyHMMER is required for this test")
    def test_from_pyhmmer(self):
        TextSequence = pyhmmer.easel.TextSequence
        sequences = [
            TextSequence(sequence="-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII", name=b"Sp8"),
            TextSequence(sequence="-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII", name=b"Sp10"),
            TextSequence(sequence="AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI", name=b"Sp26"),
            TextSequence(sequence="-----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII", name=b"Sp6"),
            TextSequence(sequence="--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF", name=b"Sp17"),
            TextSequence(sequence="-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL", name=b"Sp33"),
        ]
        alignment = Alignment.from_pyhmmer(pyhmmer.easel.TextMSA(sequences=sequences))
        self.assertEqual(list(alignment.sequences), list(self.alignment.sequences))
        self.assertEqual(list(alignment.names), list(self.alignment.names))

    @unittest.skipUnless(pyhmmer, "PyHMMER is required for this test")
    def test_to_pyhmmer(self):
        msa = self.alignment.to_pyhmmer()
        self.assertEqual(list(msa.alignment), list(self.alignment.sequences))
        self.assertEqual(list(msa.names), list(self.alignment.names))


class TestTrimmedAlignment(TestAlignment):
    def setUp(self):
        super().setUp()
        residues_mask = [True] * 46
        residues_mask[:5] = [False] * 5
        residues_mask[26:28] = [False] * 2
        sequences_mask = [True, True, False, True, True, True]
        self.trimmed = TrimmedAlignment(
            names=self.alignment.names,
            sequences=self.alignment.sequences,
            sequences_mask=sequences_mask,
            residues_mask=residues_mask,
        )

    def test_original_alignment(self):
        original = self.trimmed.original_alignment()
        self.assertEqual(original.names, self.alignment.names)
        self.assertEqual(list(original.sequences), list(self.alignment.sequences))

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
        self.assertEqual(
            self.trimmed.sequences[3], "TAPDLL-LIGFLLKTVA-TFGDTWFQLWQGLDLNKMPVF"
        )
        self.assertEqual(
            self.trimmed.sequences[-1], "--PTILNIAGLHMETDI-NFSLAWFQAWGGLEINKQAIL"
        )
        with self.assertRaises(IndexError):
            self.trimmed.sequences[5]
        with self.assertRaises(IndexError):
            self.trimmed.sequences[100]
        with self.assertRaises(IndexError):
            self.trimmed.sequences[-100]

    def test_residues_mask(self):
        mask = self.trimmed.residues_mask
        original = self.trimmed.original_alignment()
        self.assertEqual(len(mask), len(original.residues))
        self.assertEqual(
            self.trimmed.sequences[0],
            "".join([x for x, c in zip(original.sequences[0], mask) if c]),
        )

    def test_sequences_mask(self):
        mask = self.trimmed.sequences_mask
        original = self.trimmed.original_alignment()
        self.assertEqual(len(mask), len(original.sequences))
