from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from trimal cimport SequenceTypes
from trimal.alignment cimport Alignment


cdef extern from "FormatHandling/FormatManager.h" namespace "FormatHandling" nogil:

    cdef cppclass FormatManager:
        FormatManager()
        Alignment* loadAlignment(const string& inFile) except? NULL
