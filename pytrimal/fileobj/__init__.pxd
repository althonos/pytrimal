from iostream cimport streambuf

cdef extern from "pyfilebuf.h" nogil:

    cdef cppclass pyfilebuf(streambuf):
        pyfilebuf(object)
