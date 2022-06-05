from trimal.alignment cimport Alignment
from trimal.statistics cimport Similarity



cdef extern from "impl/sse.h" namespace "statistics" nogil:

    cdef cppclass SSESimilarity(Similarity):

         SSESimilarity(Alignment * parentAlignment)
