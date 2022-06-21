from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from iostream cimport istream, ostream

from trimal cimport SequenceTypes
from trimal.alignment cimport Alignment


cdef extern from "FormatHandling/FormatManager.h" namespace "FormatHandling" nogil:

    cdef cppclass FormatManager:
        FormatManager()
        Alignment* loadAlignment(const string& inFile) except? NULL
        Alignment* loadAlignment(istream& input) except? NULL
        BaseFormatHandler* getFormatFromFile(const string& filename) except? NULL
        BaseFormatHandler* getFormatFromToken(const string& token)


cdef extern from "FormatHandling/BaseFormatHandler.h" namespace "FormatHandling" nogil:

    cdef cppclass BaseFormatHandler:
        bool canLoad
        bool canSave
        string name
        string extension

        int CheckAlignment(istream *origin)
        bool SaveAlignment(const Alignment& alignment, ostream* output) except? False
        Alignment* LoadAlignment(const string& inFile) except? NULL
        Alignment* LoadAlignment(istream& input) except? NULL
