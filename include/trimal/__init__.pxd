from libcpp.string cimport string

cdef extern from "defines.h" nogil:

    cdef enum SequenceTypes:
        NotDefined
        DNA
        RNA
        AA
        DEG

cdef extern from "residueValues.h" nogil:

    const string nucleotideResidues
    const string degenerateNucleotideResidues
    const string aminoAcidResidues