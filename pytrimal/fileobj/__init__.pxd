from iostream cimport streambuf

cdef extern from "pyfilebuf.h" nogil:

    cdef cppclass pywritebuf(streambuf):
        pywritebuf(object)

    cdef cppclass pyreadbuf(streambuf):
        pyreadbuf(object)
