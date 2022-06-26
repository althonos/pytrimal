from trimal.alignment cimport Alignment
from trimal.statistics cimport Similarity


cdef extern from "impl/generic.h" namespace "statistics" nogil:
    cdef cppclass GenericSimilarity(Similarity):
         GenericSimilarity(Alignment * parentAlignment)
