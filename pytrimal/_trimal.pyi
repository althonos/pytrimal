# --- Python imports ---------------------------------------------------------

import os
import typing
from typing import BinaryIO, Sequence, List, Optional, Union, Sequence

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

# --- Constants --------------------------------------------------------------

TRIMMER_BACKEND = Literal["detect", "sse", None]
AUTOMATIC_TRIMMER_METHODS = Literal["strict", "strictplus", "gappyout", "nogaps", "noallgaps", "automated1"]

FORMATS_LOAD = Literal["clustal", "fasta", "nexus", "phylip", "phylip32", "phylip40", "pir"]
FORMATS_DUMP = Literal["clustal", "fasta", "html", "mega", "nexus", "phylip", "phylip32", "phylip40", "phylippaml", "nbrf", "pir", "fasta_m10", "nexus_m10", "phylippaml_m10", "phylip32_m10", "phylip_m10", "phylip40_m10"]

# --- Alignment classes ------------------------------------------------------

class AlignmentSequences:
    def __init__(self, alignment: Alignment) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...


class AlignmentResidues:
    def __init__(self, alignment: Alignment) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> str: ...


class Alignment:
    @typing.overload
    @classmethod
    def load(cls, file: Union[str, bytes, os.PathLike[str]], format: Optional[FORMATS_LOAD] = None) -> Alignment: ...
    @typing.overload
    @classmethod
    def load(cls, file: BinaryIO, format: FORMATS_LOAD) -> Alignment: ...
    def dump(self, file: Union[str, bytes, os.PathLike[str], BinaryIO], format: FORMATS_DUMP = "fasta") -> None: ...
    def dumps(self, format: str = "fasta", encoding: str = "utf-8") -> str: ...
    def __init__(self, names: Sequence[bytes], sequences: Sequence[str]) -> None: ...
    def __repr__(self) -> str: ...
    def __copy__(self) -> Alignment: ...
    @property
    def names(self) -> List[bytes]: ...
    @property
    def sequences(self) -> AlignmentSequences: ...
    @property
    def residues(self) -> AlignmentResidues: ...
    def copy(self) -> Alignment: ...


class TrimmedAlignment(Alignment):
    @typing.overload
    @classmethod
    def load(cls, file: Union[str, bytes, os.PathLike[str]], format: Optional[FORMATS_LOAD] = None) -> TrimmedAlignment: ...
    @typing.overload
    @classmethod
    def load(cls, file: BinaryIO, format: FORMATS_LOAD) -> TrimmedAlignment: ...
    def __init__(
        self,
        names: Sequence[bytes],
        sequences: Sequence[str],
        sequences_mask: Optional[Sequence[bool]] = None,
        residues_mask: Optional[Sequence[bool]] = None,
    ) -> None: ...
    def original_alignment(self) -> Alignment: ...
    def terminal_only(self) -> TrimmedAlignment: ...
    def copy(self) -> TrimmedAlignment: ...
    @property
    def residues_mask(self) -> List[bool]: ...
    @property
    def sequences_mask(self) -> List[bool]: ...


# -- Trimmer classes ---------------------------------------------------------

class BaseTrimmer:
    def __init__(self, *, backend: TRIMMER_BACKEND = "detect") -> None: ...
    def trim(self, alignment: Alignment, matrix: Optional[SimilarityMatrix] = None) -> TrimmedAlignment: ...


class AutomaticTrimmer(BaseTrimmer):
    def __init__(self, method: AUTOMATIC_TRIMMER_METHODS = "strict", *, backend: TRIMMER_BACKEND = "detect") -> None: ...


class ManualTrimmer(BaseTrimmer):
    def __init__(
        self,
        *,
        gap_threshold: Optional[float] = None,
        gap_absolute_threshold: Optional[int] = None,
        similarity_threshold: Optional[float] = None,
        consistency_threshold: Optional[float] = None,
        conservation_percentage: Optional[float] = None,
        window: Optional[int] = None,
        gap_window: Optional[int] = None,
        similarity_window: Optional[int] = None,
        consistency_window: Optional[int] = None,
        backend: TRIMMER_BACKEND = "detect",
    ) -> None: ...


# -- Misc classes ------------------------------------------------------------


class SimilarityMatrix:
    @classmethod
    def aa(cls) -> SimilarityMatrix: ...
    @classmethod
    def nt(cls, degenerated: bool =False) -> SimilarityMatrix: ...
    def __init__(self, alphabet: str, matrix: Sequence[Sequence[float]]) -> None: ...
    def __len__(self) -> int: ...
    def distance(self, a: str, b: str) -> float: ...
    def similarity(self, a: str, b: str) -> float: ...
