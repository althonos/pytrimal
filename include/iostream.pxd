from libcpp.string cimport string

cdef extern from "<streambuf>" namespace "std" nogil:
    cdef cppclass streambuf:
        pass

cdef extern from "<fstream>" namespace "std" nogil:
    cdef cppclass filebuf(streambuf):
        pass

cdef extern from "<sstream>" namespace "std" nogil:
    cdef cppclass stringbuf(streambuf):
        stringbuf()
        stringbuf(const string& str)

        string str() const
        void str (const string& str)

cdef extern from "<ostream>" namespace "std" nogil:
    cdef cppclass ostream:
        ostream(streambuf* sb)
