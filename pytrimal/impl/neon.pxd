from trimal.alignment cimport Alignment
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Similarity, Gaps


cdef extern from "impl/neon.h" namespace "statistics" nogil:
    cdef cppclass NEONSimilarity(Similarity):
         NEONSimilarity(Alignment * parentAlignment)
    cdef cppclass NEONGaps(Gaps):
         NEONGaps(Alignment* parentAlignment)

cdef extern from "impl/sse.h" nogil:
    cdef cppclass NEONCleaner(Cleaner):
        NEONCleaner(Alignment* parentAlignment)
