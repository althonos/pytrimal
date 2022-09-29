from libcpp cimport bool

from trimal.similarity_matrix cimport similarityMatrix


cdef extern from "Alignment/Alignment.h" nogil:

    cdef cppclass Alignment:
        pass


cdef extern from "Statistics/Consistency.h" namespace "statistics" nogil:

    cdef cppclass Consistency:
        pass


cdef extern from "Statistics/Gaps.h" namespace "statistics" nogil:

    cdef cppclass Gaps:
        void CalculateVectors()


cdef extern from "Statistics/Mold.h" namespace "statistics" nogil:

    cdef cppclass Mold:
        pass


cdef extern from "Statistics/Similarity.h" namespace "statistics" nogil:

    cdef cppclass Similarity:
        Alignment* alig
        int halfWindow
        float* MDK
        float* MDK_Window
        float **matrixIdentity
        similarityMatrix* simMatrix
        int* refCounter

        Similarity(Alignment * parentAlignment)
        Similarity(Alignment* parentAlignment, Similarity* mold)
        bool calculateVectors(bool cutByGap = true)
        void calculateMatrixIdentity()
        bool applyWindow(int halfW)
        bool isDefinedWindow()
        float *getMdkWindowedVector()
        bool setSimilarityMatrix(similarityMatrix * sm)
        bool isSimMatrixDef()
        double calcCutPoint(float baseLine, float conservationPct)
        void printConservationColumns()
        void printConservationAcl()


cdef extern from "Statistics/Manager.h" namespace "statistics" nogil:

    cdef cppclass Manager:
        Gaps *gaps
        Similarity *similarity
        Consistency *consistency
        similarityMatrix* _similarityMatrix
        int ghWindow
        int shWindow
        Alignment *alig

        Manager(Alignment *parent)
        Manager(Alignment *parent, Manager *mold)
        bool setSimilarityMatrix(similarityMatrix *sm)
        bool calculateGapStats() except? False
        void printStatisticsGapsColumns()
        void printStatisticsGapsTotal()
        bool calculateConservationStats() except? False
        void printStatisticsConservationColumns()
        void printStatisticsConservationTotal()
        void printCorrespondence()
