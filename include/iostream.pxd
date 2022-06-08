from libcpp.string cimport string

cdef extern from "<ios>" namespace "std" nogil:

    cdef cppclass ios_base:
        cppclass openmode:
            openmode operator|(openmode)

cdef extern from "<streambuf>" namespace "std" nogil:
    cdef cppclass streambuf:
        char* eback()
        char* egptr()

cdef extern from "<fstream>" namespace "std" nogil:
    cdef cppclass filebuf(streambuf):
        filebuf()
        filebuf* open(const char* filename,  ios_base.openmode mode);

cdef extern from "<sstream>" namespace "std" nogil:
    cdef cppclass stringbuf(streambuf):
        stringbuf()
        stringbuf(const string& str)

        string str() const
        void str (const string& str)

cdef extern from "<ostream>" namespace "std" nogil:
    cdef cppclass ostream:
        ostream(streambuf* sb)
