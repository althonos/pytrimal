from trimal.alignment cimport Alignment
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Similarity, Gaps


cdef extern from "impl/generic.h" namespace "statistics" nogil:
    cdef cppclass GenericSimilarity(Similarity):
        GenericSimilarity(Alignment * parentAlignment)
    cdef cppclass GenericGaps(Gaps):
        GenericGaps(Alignment* parentAlignment)


cdef extern from "impl/sse.h" nogil:
    cdef cppclass GenericCleaner(Cleaner):
        GenericCleaner(Alignment* parentAlignment)
