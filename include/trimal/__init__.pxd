cdef extern from "defines.h" nogil:

    cdef enum SequenceTypes:
        NotDefined
        DNA
        RNA
        AA
        DEG
