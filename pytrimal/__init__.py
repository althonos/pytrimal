# noqa: D104
from ._version import __version__  # isort: skip

from . import _trimal
from ._trimal import (
    Alignment,
    AlignmentSequences,
    AutomaticTrimmer,
    BaseAlignment,
    BaseTrimmer,
    ManualTrimmer,
    SimilarityMatrix,
    TrimmedAlignment,
)

__doc__ = _trimal.__doc__
__all__ = [
    "Alignment",
    "AlignmentSequences",
    "BaseAlignment",
    "TrimmedAlignment",
    "BaseTrimmer",
    "AutomaticTrimmer",
    "ManualTrimmer",
    "SimilarityMatrix"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
