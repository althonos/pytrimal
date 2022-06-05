from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from trimal cimport SequenceTypes
from trimal.alignment cimport Alignment


cdef extern from "trimalManager.h" nogil:

    cdef cppclass trimAlManager:
        vector[string]* vcfs

        bool appearErrors
        bool getComplementary
        bool getComplementarySeq
        bool columnNumbering
        bool nogaps
        bool noallgaps
        bool gappyout
        bool strict
        bool strictplus
        bool automated1
        bool sgc
        bool sgt
        bool ssc
        bool sst
        bool sfc
        bool sft
        bool sident
        bool soverlap
        bool selectSeqs
        bool selectCols
        bool splitByStopCodon
        bool terminalOnly
        bool keepSeqs
        bool ignoreStopCodon
        bool ignoreFilter
        bool removeDuplicates

        float conservationThreshold
        float gapThreshold
        float similarityThreshold
        float consistencyThreshold
        float residuesOverlap
        float sequenceOverlap
        float maxIdentity
        float minCoverage
        float minQuality

        int i
        int stats
        int windowSize
        int gapWindow
        int similarityWindow
        int consistencyWindow
        int blockSize
        int clusters
        int automatedMethodCount
        int alternative_matrix
        int gapAbsoluteThreshold

        int* delColumns
        int* delSequences
        int* sequencesLengths

        string* sequencesNames

        Alignment* origAlig
        Alignment* singleAlig
        Alignment* tempAlig
        Alignment* backtranslationAlig
        Alignment** compareAlignmentsArray

        trimAlManager()

        int perform()
        int parseArguments(int argc, char** argv) except *
        int processArguments(char** argv) except *

        bool performCompareset()
        void output_reports()
        void save_alignment()
        void svg_stats_out()
        void print_statistics()
        bool create_or_use_similarity_matrix() except False
        void clean_alignment() except *
        void postprocess_alignment()
        void set_window_size()
