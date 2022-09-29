from trimal.alignment cimport Alignment
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Similarity, Gaps


cdef extern from "impl/sse.h" namespace "statistics" nogil:
    cdef cppclass SSESimilarity(Similarity):
        SSESimilarity(Alignment* parentAlignment)
    cdef cppclass SSEGaps(Gaps):
        SSEGaps(Alignment* parentAlignment)


cdef extern from "impl/sse.h" nogil:
    cdef cppclass SSECleaner(Cleaner):
        SSECleaner(Alignment* parentAlignment)
