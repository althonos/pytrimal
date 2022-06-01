# distutils: language = c++
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libc.math cimport isnan, NAN
from libcpp cimport bool
from libcpp.string cimport string

cimport trimal
cimport trimal.alignment
cimport trimal.format_manager
cimport trimal.manager
cimport trimal.report_system
cimport trimal.similarity_matrix


# --- Python imports ---------------------------------------------------------

import os
import threading


# --- Constants --------------------------------------------------------------

cdef set AUTOMATED_TRIMMER_METHODS = {
    "strict",
    "strictplus",
    "gappyout",
    "nogaps",
    "noallgaps",
    "automated1",
}


# --- Utilities --------------------------------------------------------------

cdef float _check_range(object value, str name, float min_value, float max_value) except NAN:
    if value < min_value or value > max_value or isnan(value):
        raise ValueError(f"Invalid value for `{name}`: {value!r}")
    return value


# --- Cython classes ---------------------------------------------------------


cdef class Alignment:

    cdef trimal.alignment.Alignment* _ali

    def __cinit__(self):
        self._ali = NULL

    def __dealloc__(self):
        if self._ali != NULL:
            self._ali.saveResidues = NULL
            # del self._ali

    @classmethod
    def load(cls, object path not None):
        cdef trimal.format_manager.FormatManager manager
        cdef Alignment alignment = cls.__new__(cls)
        cdef string    path_ =  os.fsencode(path)
        alignment._ali = manager.loadAlignment(path_)
        return alignment

    @property
    def names(self):
        return [
            self._ali.seqsName[i]
            for i in range(self._ali.originalNumberOfSequences)
            if self._ali.saveSequences[i] != -1
        ]

    @property
    def sequences(self):
        return [
            bytes([
                self._ali.sequences[i][j]
                for j in range(self._ali.originalNumberOfResidues)
                if self._ali.saveResidues is NULL or self._ali.saveResidues[j] != -1
            ]).decode()
            for i in range(self._ali.originalNumberOfSequences)
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1
        ]


cdef class TrimmedAlignment(Alignment):
    pass


cdef class BaseTrimmer:

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        pass

    cpdef TrimmedAlignment trim(self, Alignment alignment):
        # use a local manager object so that this method is re-entrant
        cdef trimal.manager.trimAlManager _mg

        # copy the alignment to the manager object so that the original
        # alignment is left untouched
        _mg.origAlig = new trimal.alignment.Alignment(alignment._ali[0])

        # configure the manager (to be implemented by the different subclasses)
        self._configure_manager(&_mg)

        # set flags
        # self._mg.origAlig.Cleaning.setTrimTerminalGapsFlag(self.terminal_only)
        # self._mg.origAlig.setKeepSequencesFlag(self.keep_sequences)
        _mg.set_window_size()
        if _mg.blockSize != -1:
            _mg.origAlig.setBlockSize(self._mg.blockSize)

        _mg.create_or_use_similarity_matrix()
        # self._mg.print_statistics()
        _mg.clean_alignment()

        if _mg.singleAlig == NULL:
            _mg.singleAlig = _mg.origAlig
            _mg.origAlig = NULL

        _mg.postprocess_alignment()
        # _mg.output_reports()
        # _mg.save_alignment()

        cdef TrimmedAlignment trimmed = TrimmedAlignment.__new__(TrimmedAlignment)
        trimmed._ali = new trimal.alignment.Alignment(_mg.singleAlig[0])
        return trimmed


cdef class AutomaticTrimmer(BaseTrimmer):
    """A sequence alignment trimmer with automatic parameter detection.
    """

    cdef readonly str method

    def __init__(self, str method="strict"):
        """__init__(self, method="strict")\n--

        Create a new automatic alignment trimmer using the given method.

        Arguments:
            method (`str`): The automatic aligment trimming method.

        Raises:
            `ValueError`: When ``method`` is not one of the automatic
                alignment trimming methods supported by trimAl.

        """
        super().__init__()

        if method not in AUTOMATED_TRIMMER_METHODS:
            raise ValueError(f"Invalid value for `method`: {method!r}")
        self.method = method

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount = 1
        if self.method == "strict":
            manager.strict = True
        elif self.method == "strictplus":
            manager.strictplus = True
        elif self.method == "gappyout":
            manager.gappyout = True
        elif self.method == "nogaps":
            manager.nogaps = True
        elif self.method == "noallgaps":
            manager.noallgaps = True
        elif self.method == "automated1":
            manager.automated1 = True


cdef class ManualTrimmer(BaseTrimmer):
    """A sequence alignment trimmer with manually defined thresholds.
    """

    cdef float   _gap_threshold
    cdef ssize_t _gap_absolute_threshold
    cdef float   _similarity_threshold
    cdef float   _consistency_threshold
    cdef float   _conservation_percentage

    def __cinit__(self):
        self._gap_threshold           = -1
        self._gap_absolute_threshold  = -1
        self._similarity_threshold    = -1
        self._consistency_threshold   = -1
        self._conservation_percentage = -1

    def __init__(
        self,
        *,
        object gap_threshold           = None,
        object gap_absolute_threshold  = None,
        object similarity_threshold    = None,
        object consistency_threshold   = None,
        object conservation_percentage = None,
    ):
        """__init__(self, *, gap_threshold=None, gap_absolute_threshold=None, similarity_threshold=None, consistency_threshold=None, conservation_percentage=None)\n--

         Create a new manual alignment trimmer with the given parameters.

         Keyword Arguments:
             gap_threshold (`float`, *optional*): The minimum fraction of
                 non-gap characters that must be present in a column to
                 keep the column.
            gap_absolute_threshold (`int`, *optional*): The absolute number
                of gaps allowed on a column to keep it in the alignment.
                Incompatible with ``gap_threshold``.
            similarity_threshold (`float`, *optional*): The minimum average
                similarity required.
            consistency_threshold (`float`, *optional*): The minimum
                consistency value required.
            conservation_percentage (`float`, *optional*): The minimum
                percentage of positions in the original alignment to
                conserve.

        """
        super().__init__()

        if gap_threshold is not None and gap_absolute_threshold is not None:
            raise ValueError("Cannot specify both `gap_threshold` and `gap_absolute_threshold`")

        if gap_threshold is not None:
            self._gap_threshold = 1 - _check_range(gap_threshold, "gap_threshold", 0, 1)
        if gap_absolute_threshold is not None:
            if gap_absolute_threshold < 0 or gap_absolute_threshold > 100:
                raise ValueError(f"Invalid value for `gap_absolute_threshold`: {gap_absolute_threshold!r}")
            self._gap_absolute_threshold = gap_absolute_threshold
        if similarity_threshold is not None:
            self._similarity_threshold = _check_range(similarity_threshold, "similarity_threshold", 0, 1)
        if consistency_threshold is not None:
            self._consistency_threshold = _check_range(consistency_threshold, "consistency_threshold", 0, 1)
        if conservation_percentage is not None:
            self._conservation_percentage = _check_range(conservation_percentage, "conservation_percentage", 0, 100)

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount = 0
        manager.gapThreshold = self._gap_threshold
        manager.gapAbsoluteThreshold = self._gap_absolute_threshold
        manager.similarityThreshold = self._similarity_threshold
        manager.consistencyThreshold = self._consistency_threshold
        manager.conservationThreshold = self._conservation_percentage


cdef class SimilarityMatrix:

    cdef trimal.similarity_matrix.similarityMatrix _smx

    @classmethod
    def aa(cls):
        cdef SimilarityMatrix matrix = cls.__new__(cls)
        with nogil:
            matrix._smx.defaultAASimMatrix()
        return matrix

    @classmethod
    def nt(cls, bool degenerated=False):
        cdef SimilarityMatrix matrix = cls.__new__(cls)
        with nogil:
            if degenerated:
                matrix._smx.defaultNTSimMatrix()
            else:
                matrix._smx.defaultNTDegeneratedSimMatrix()
        return matrix

    def __cinit__(self):
        self._smx

    def __init__(self):
        raise NotImplementedError("SimilarityMatrix.__init__")

    cpdef float distance(self, str a, str b) except -1:
        """distance(self, a, b)\n--

        Return the distance between two sequence characters.

        """
        cdef float distance = 0.0
        cdef char  a_code   = ord(a)
        cdef char  b_code   = ord(b)
        with nogil:
            distance  = self._smx.getDistance(a_code, b_code)
        return distance
