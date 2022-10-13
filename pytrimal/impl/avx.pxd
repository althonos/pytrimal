from trimal.alignment cimport Alignment
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Similarity, Gaps


cdef extern from "impl/avx.h" namespace "statistics" nogil:
    cdef cppclass AVXSimilarity(Similarity):
         AVXSimilarity(Alignment * parentAlignment)
    cdef cppclass AVXGaps(Gaps):
         AVXGaps(Alignment* parentAlignment)

cdef extern from "impl/avx.h" nogil:
    cdef cppclass AVXCleaner(Cleaner):
        AVXCleaner(Alignment* parentAlignment)
