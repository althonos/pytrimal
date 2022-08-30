# noqa: D104
from ._version import __version__  # isort: skip

from . import _trimal
from ._trimal import (
    Alignment,
    AlignmentResidues,
    AlignmentSequences,
    AutomaticTrimmer,
    BaseTrimmer,
    ManualTrimmer,
    OverlapTrimmer,
    RepresentativeTrimmer,
    SimilarityMatrix,
    TrimmedAlignment,
)

__doc__ = _trimal.__doc__
__all__ = [
    "Alignment",
    "AlignmentResidues",
    "AlignmentSequences",
    "TrimmedAlignment",
    "BaseTrimmer",
    "AutomaticTrimmer",
    "ManualTrimmer",
    "OverlapTrimmer",
    "SimilarityMatrix",
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pytrimal.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
