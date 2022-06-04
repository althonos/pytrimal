from libcpp cimport bool
from libcpp.string cimport string

from trimal cimport SequenceTypes
from trimal.cleaner cimport Cleaner
from trimal.statistics cimport Manager


cdef extern from "Alignment/Alignment.h" nogil:

    cdef cppclass Alignment:
        int dataType

        Cleaner* Cleaning
        Manager* Statistics
        int* SeqRef
        int originalNumberOfSequences
        int numberOfSequences
        int originalNumberOfResidues
        int numberOfResidues
        bool isAligned
        string* sequences
        string* seqsName
        string* seqsInfo
        string filename
        string alignmentInfo
        float** identities
        float** overlaps
        int* saveResidues
        int* saveSequences

        cppclass sequencesMatrix:
            int resNumber
            int seqsNumber
            int **matrix
            string *seqsName

        bool fillMatrices(bool aligned, bool checkInvalidChars) except False

        Alignment()
        Alignment(Alignment&)

        int getNumSpecies()
        void getSequences(string* names)
        void getSequences(string* names, int* lenghts)
        void getSequences(string* names, string* sequences, int* lenghts)
        bool getSequenceNameOrder(string* names, int* orderVector)
        int getNumAminos()

        void setWindowsSize(int ghWindow, int shWindow)
        void setBlockSize(int blockSize)
        void calculateSeqOverlap()
        void printSeqIdentity()
        void printSeqOverlap()
        int getAlignmentType() const

        bool isFileAligned()
        Alignment *getTranslationCDS(Alignment *proteinAlignment)
        bool checkCorrespondence(string *names, int *lenghts, int totalInputSequences, int multiple)
        void calculateColIdentity(float *columnIdentity)
        void setKeepSequencesFlag(bool newFlagValue)

        # void printAlignmentInfo(ostream &output)

        bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, Alignment *proteinAlignment)
        bool alignmentSummaryHTML(const Alignment &trimmedAlig, const char * const destFile)
        bool statSVG(const char *const destFile)
        bool alignmentSummarySVG(Alignment &trimmedAlig, const char *destFile, int blocks)
        void updateSequencesAndResiduesNums(bool countSequences = true, bool countResidues = true)
