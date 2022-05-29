# distutils: language = c++
# cython: language_level=3, linetrace=True

from libcpp cimport bool

cimport trimal
cimport trimal.report_system
from trimal.similarity_matrix cimport similarityMatrix

# --- Error management -------------------------------------------------------


cdef class SimilarityMatrix:

    cdef similarityMatrix _smx

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
