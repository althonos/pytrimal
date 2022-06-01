from ._version import __version__

from . import _trimal
from ._trimal import (
    Alignment,
    TrimmedAlignment,
    BaseTrimmer,
    AutomaticTrimmer,
    ManualTrimmer,
    SimilarityMatrix
)

__doc__ = _trimal.__doc__
__all__ = [
    "Alignment",
    "TrimmedAlignment",
    "BaseTrimmer",
    "AutomaticTrimmer",
    "ManualTrimmer",
    "SimilarityMatrix"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
