cdef extern from "Statistics/similarityMatrix.h" namespace "statistics":

    cdef cppclass similarityMatrix:
        int* vhash
        float** simMat
        float** distMat
        int numPositions

        void defaultAASimMatrix() nogil except +
        void defaultNTSimMatrix() nogil except +
        void defaultNTDegeneratedSimMatrix() nogil except +

        float getDistance(char a, char b) except -1
