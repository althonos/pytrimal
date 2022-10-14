from trimal.alignment cimport Alignment
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Similarity, Gaps


cdef extern from "impl/mmx.h" namespace "statistics" nogil:
    cdef cppclass MMXSimilarity(Similarity):
        MMXSimilarity(Alignment* parentAlignment)
    cdef cppclass MMXGaps(Gaps):
        MMXGaps(Alignment* parentAlignment)


cdef extern from "impl/mmx.h" nogil:
    cdef cppclass MMXCleaner(Cleaner):
        MMXCleaner(Alignment* parentAlignment)
