from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

from trimal cimport SequenceTypes
from trimal.alignment cimport Alignment

cdef extern from "FormatHandling/FormatManager.h" namespace "FormatHandling" nogil:

    cdef cppclass FormatManager:
        FormatManager()
        Alignment* loadAlignment(const string& inFile) except? NULL
