from iostream cimport streambuf

cdef extern from "pywritebuf.h" nogil:

    cdef cppclass pywritebuf(streambuf):
        pywritebuf(object)


cdef extern from "pyreadbuf.h" nogil:

    cdef cppclass pyreadbuf(streambuf):
        pyreadbuf(object)


cdef extern from "pyreadintobuf.h" nogil:

    cdef cppclass pyreadintobuf(streambuf):
        pyreadintobuf(object)
