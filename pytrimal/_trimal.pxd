# distutils: language = c++
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libcpp cimport bool

cimport trimal
cimport trimal.alignment
cimport trimal.manager
cimport trimal.similarity_matrix


# --- Alignment classes ------------------------------------------------------

cdef class AlignmentSequences:
    cdef trimal.alignment.Alignment* _ali            # pointer to raw alignment
    cdef Alignment                   _owner          # reference to owner
    cdef int*                        _index_mapping  # old-to-new index mapping, or NULL
    cdef ssize_t                     _length         # number of sequences
    cdef bool                        _free_mapping   # whether to free the _index_mapping

    cdef str _sequence(self, int index)
    cdef AlignmentSequences _slice(self, int start, int stop, int stride)


cdef class AlignmentResidues:
    cdef trimal.alignment.Alignment* _ali           # pointer to raw alignment
    cdef Alignment                   _owner         # reference to owner
    cdef int*                        _index_mapping # old-to-new index mapping, or NULL
    cdef ssize_t                     _length        # number of residues
    cdef bool                        _free_mapping  # whether to free the _index_mapping

    cdef str _column(self, int index)
    cdef AlignmentResidues _slice(self, int start, int stop, int stride)


cdef class Alignment:
    cdef trimal.alignment.Alignment* _ali
    cdef int*                        _sequences_mapping
    cdef int*                        _residues_mapping

    cpdef Alignment copy(self)
    cpdef str dumps(self, str format=*, str encoding=*)
    cpdef void dump(self, object file, str format=*) except *


cdef class TrimmedAlignment(Alignment):
    cdef void _build_index_mapping(self) except *
    cpdef Alignment original_alignment(self)
    cpdef TrimmedAlignment terminal_only(self)
    cpdef TrimmedAlignment copy(self)


# -- Trimmer classes ---------------------------------------------------------

cdef class BaseTrimmer:
    cdef int _backend

    cdef void _setup_simd_code(self, trimal.manager.trimAlManager* manager) nogil
    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)
    cpdef TrimmedAlignment trim(self, Alignment alignment, SimilarityMatrix matrix = ?)


cdef class AutomaticTrimmer(BaseTrimmer):
    cdef readonly str method

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)


cdef class ManualTrimmer(BaseTrimmer):
    cdef float   _gap_threshold
    cdef ssize_t _gap_absolute_threshold
    cdef float   _similarity_threshold
    cdef float   _consistency_threshold
    cdef float   _conservation_percentage
    cdef int     _window
    cdef int     _gap_window
    cdef int     _similarity_window
    cdef int     _consistency_window

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)


cdef class OverlapTrimmer(BaseTrimmer):
    cdef float _sequence_overlap
    cdef float _residue_overlap

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)


cdef class RepresentativeTrimmer(BaseTrimmer):
    cdef int    _clusters
    cdef float  _identity_threshold

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)


# -- Misc classes ------------------------------------------------------------


cdef class SimilarityMatrix:
    cdef Py_ssize_t                                _suboffsets[2]
    cdef Py_ssize_t                                _shape[2]
    cdef Py_ssize_t                                _strides[2]
    cdef trimal.similarity_matrix.similarityMatrix _smx

    cpdef float similarity(self, str a, str b) except -1
    cpdef float distance(self, str a, str b) except -1
