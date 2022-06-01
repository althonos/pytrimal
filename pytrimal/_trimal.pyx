# distutils: language = c++
# cython: language_level=3, linetrace=True
"""Bindings to trimAl, a tool for automated alignment trimming.

References:
    - Capella-Gutiérrez, Salvador, José M. Silla-Martínez, and Toni Gabaldón.
      *TrimAl: A Tool for Automated Alignment Trimming in Large-Scale
      Phylogenetic Analyses*. Bioinformatics 25, no. 15 (2009): 1972–73.
      :doi:`10.1093/bioinformatics/btp348`.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.ref cimport Py_INCREF

from _unicode cimport (
    Py_UCS1,
    PyUnicode_1BYTE_DATA,
    PyUnicode_1BYTE_KIND,
    PyUnicode_New,
    PyUnicode_READY,
)
from libc.math cimport NAN, isnan
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

@cython.freelist(32)
cdef class _AlignmentSequences:
      """A read-only view over the sequences of an alignment.
      """

      cdef trimal.alignment.Alignment* _ali
      cdef Alignment                   _owner
      cdef int*                        _sequence_mapping

      def __cinit__(self, Alignment alignment):
          self._owner = alignment
          self._ali = alignment._ali
          self._sequence_mapping = alignment._sequence_mapping

      def __len__(self):
          assert self._ali is not NULL
          return self._ali.numberOfSequences

      def __getitem__(self, ssize_t index):
          assert self._ali is not NULL

          cdef ssize_t index_ = index
          cdef ssize_t length = self._ali.numberOfSequences

          if index_ < 0:
              index_ += length
          if index_ < 0 or index_ > length:
              raise IndexError(index)

          cdef size_t   i
          cdef size_t   x         = 0
          cdef int      old_index = self._sequence_mapping[index_]
          cdef str      seq       = PyUnicode_New(self._ali.numberOfResidues, 0x7f)

          IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
              PyUnicode_READY(seq)

          cdef Py_UCS1* seqdata = PyUnicode_1BYTE_DATA(seq)

          for i in range(self._ali.originalNumberOfResidues):
              if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                  seqdata[x] = self._ali.sequences[old_index][i]
                  x += 1

          return seq

          # return _AlignmentSequence.__new__(_AlignmentSequence, self._owner, index_)


cdef class Alignment:
    """A multiple sequence alignment.
    """

    cdef trimal.alignment.Alignment* _ali
    cdef int*                        _sequence_mapping

    def __cinit__(self):
        self._ali = NULL

    def __dealloc__(self):
        if self._sequence_mapping != NULL:
            PyMem_Free(self._sequence_mapping)

    @classmethod
    def load(cls, object path not None):
        cdef trimal.format_manager.FormatManager manager
        cdef Alignment alignment = cls.__new__(cls)
        cdef string    path_ =  os.fsencode(path)
        alignment._ali = manager.loadAlignment(path_)
        return alignment

    @property
    def names(self):
        """sequence of `bytes`: The names of the sequences in the alignment.
        """
        assert self._ali is not NULL

        cdef size_t i
        cdef size_t x     = 0
        cdef bytes  name
        cdef object names = PyList_New(self._ali.numberOfSequences)

        for i in range(self._ali.originalNumberOfSequences):
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                name = PyBytes_FromStringAndSize( self._ali.seqsName[i].data(), self._ali.seqsName[i].size() )
                PyList_SET_ITEM(names, x, name)
                Py_INCREF(name)  # manually increase reference count because `PyList_SET_ITEM` doesn't
                x += 1

        return names

    @property
    def sequences(self):
        """sequence of `str`: The sequences in the alignment.
        """
        assert self._ali is not NULL

        cdef int i
        cdef int x = 0

        # build a mapping of old index to new index for sequences in the
        # alignment ()
        if self._sequence_mapping is NULL:
            self._sequence_mapping = <int*> PyMem_Malloc(self._ali.numberOfSequences * sizeof(int))
            for i in range(self._ali.originalNumberOfSequences):
                if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                    self._sequence_mapping[x] = i
                    x += 1

        return _AlignmentSequences(self)


cdef class TrimmedAlignment(Alignment):
    pass


cdef class BaseTrimmer:
    """A sequence alignment trimmer.

    All subclasses provide the same `trim` method, and are configured
    through their constructor.

    """

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        pass

    cpdef TrimmedAlignment trim(self, Alignment alignment):
        """trim(self, alignment)\n--

        Trim the provided alignment.

        Returns:
            `~pytrimal.TrimmedAlignment`: The trimmed alignment.

        Note:
            This method is re-entrant, and can be called safely accross
            different threads.

        """
        # use a local manager object so that this method is re-entrant
        cdef trimal.manager.trimAlManager _mg

        # copy the alignment to the manager object so that the original
        # alignment is left untouched
        _mg.origAlig = new trimal.alignment.Alignment(alignment._ali[0])

        # configure the manager (to be implemented by the different subclasses)
        self._configure_manager(&_mg)

        with nogil:
            # set flags
            # self._mg.origAlig.Cleaning.setTrimTerminalGapsFlag(self.terminal_only)
            # self._mg.origAlig.setKeepSequencesFlag(self.keep_sequences)
            _mg.set_window_size()
            if _mg.blockSize != -1:
                _mg.origAlig.setBlockSize(_mg.blockSize)

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
        BaseTrimmer._configure_manager(self, manager)
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
        BaseTrimmer._configure_manager(self, manager)
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
