cdef extern from "Statistics/similarityMatrix.h" namespace "statistics" nogil:

    cdef cppclass similarityMatrix:
        int* vhash
        float** simMat
        float** distMat
        int numPositions

        similarityMatrix()
        void memoryAllocation(int nPos)

        void defaultAASimMatrix() except +
        void defaultNTSimMatrix() except +
        void defaultNTDegeneratedSimMatrix() except +

        float getDistance(char a, char b) except -1
