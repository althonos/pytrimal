from libcpp cimport bool
from libcpp.string cimport string

from trimal cimport SequenceTypes


cdef extern from "Alignment/sequencesMatrix.h" namespace "Alignment" nogil:

    cdef sequencesMatrix:
        int resNumber
        int seqsNumber
        int **matrix
        string* seqsName
