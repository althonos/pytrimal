try:
    try:
        from importlib.resources import files
    except ImportError:
        from importlib_resources import files
except ImportError:
    files = None

from .. import Alignment


class TrimmerTestCase(object):
    backend = None

    def assertTrimmedAlignmentEqual(self, trimmed, expected):
        self.assertEqual(len(trimmed.names), len(expected.names))
        self.assertEqual(len(trimmed.sequences), len(expected.sequences))
        self.assertEqual(trimmed.names, expected.names)
        for seq1, seq2 in zip(trimmed.sequences, expected.sequences):
            self.assertEqual(seq1, seq2)

    @classmethod
    def _open_data(cls, name, mode="r"):
        return files(__package__).joinpath("data", name).open(mode=mode)

    @classmethod
    def _load_alignment(cls, name, format="fasta"):
        with cls._open_data(name, "rb") as file:
            return Alignment.load(file, format=format)