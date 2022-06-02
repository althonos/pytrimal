from libcpp cimport bool


cdef extern from "Alignment/Alignment.h" nogil:
    cdef cppclass Alignment:
        pass

cdef extern from "Cleaner.h" nogil:
    cdef cppclass Cleaner:
        bool terminalOnly
        bool keepSequences
        int blockSize
        int left_boundary
        int right_boundary

        int selectMethod()
        Alignment *cleanByCutValueOverpass(double cut, float baseLine, const int *gInCol, bool complementary)
        Alignment *cleanByCutValueFallBehind(float cut, float baseLine, const float *ValueVect, bool complementary)
        Alignment *cleanByCutValueOverpassOrEquals(double cutGaps, const int *gInCol, float baseLine, float cutCons, const float *MDK_Win, bool complementary)
        Alignment *cleanStrict(int gapCut, const int *gInCol, float simCut, const float *MDK_W, bool complementary, bool variable)
        Alignment *cleanOverlapSeq(float minimumOverlap, float *overlapSeq, bool complementary)
        Alignment *cleanGaps(float baseLine, float gapsPct, bool complementary)
        Alignment *cleanConservation(float baseLine, float conservationPct, bool complementary)
        Alignment *clean(float baseLine, float GapsPct, float conservationPct, bool complementary)
        Alignment *cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary)

        bool calculateSpuriousVector(float overlap, float *spuriousVector)
        Alignment *cleanSpuriousSeq(float overlapColumn, float minimumOverlap, bool complementary)
        Alignment *clean2ndSlope(bool complementary)
        Alignment *cleanCombMethods(bool complementary, bool variable)
        Alignment *cleanNoAllGaps(bool complementary)
        Alignment *removeColumns(int *columns, int init, int size, bool complementary)
        Alignment *removeSequences(int *seqs, int init, int size, bool complementary)
        Alignment *getClustering(float identityThreshold)
        float getCutPointClusters(int clusterNumber)

        void removeSmallerBlocks(int blockSize, Alignment &original)
        bool removeOnlyTerminal() except False
        void removeAllGapsSeqsAndCols(bool seqs = true, bool cols = true) except *
        void setTrimTerminalGapsFlag(bool terminalOnly_)
        void setBoundaries(int *boundaries)
        void calculateSeqIdentity() except *
        void calculateRelaxedSeqIdentity()
        int *calculateRepresentativeSeq(float maximumIdent)
        void computeComplementaryAlig(bool residues, bool sequences)
        void removeDuplicates() except *
