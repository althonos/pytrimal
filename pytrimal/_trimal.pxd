# distutils: language = c++
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

cimport trimal
cimport trimal.alignment
cimport trimal.manager
cimport trimal.similarity_matrix


# --- Alignment classes ------------------------------------------------------

cdef class AlignmentSequences:
    cdef trimal.alignment.Alignment* _ali
    cdef Alignment                   _owner
    cdef int*                        _index_mapping


cdef class AlignmentResidues:
    cdef trimal.alignment.Alignment* _ali
    cdef Alignment                   _owner
    cdef int*                        _index_mapping


cdef class Alignment:
    cdef trimal.alignment.Alignment* _ali
    cdef int*                        _sequences_mapping
    cdef int*                        _residues_mapping

    cpdef Alignment copy(self)


cdef class TrimmedAlignment(Alignment):
    cdef void _build_index_mapping(self) except *
    cpdef Alignment original_alignment(self)
    cpdef TrimmedAlignment terminal_only(self)
    cpdef TrimmedAlignment copy(self)


# -- Trimmer classes ---------------------------------------------------------

cdef class BaseTrimmer:
    cdef int _backend

    cdef void _setup_simd_code(self, trimal.manager.trimAlManager* manager)
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

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager)


# -- Misc classes ------------------------------------------------------------


cdef class SimilarityMatrix:
    cdef Py_ssize_t                                _suboffsets[2]
    cdef Py_ssize_t                                _shape[2]
    cdef Py_ssize_t                                _strides[2]
    cdef trimal.similarity_matrix.similarityMatrix _smx

    cpdef float similarity(self, str a, str b) except -1
    cpdef float distance(self, str a, str b) except -1
