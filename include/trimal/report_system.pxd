from libcpp cimport bool
from libcpp.map cimport map


cdef extern from "reportsystem.h" namespace "reporting" nogil:

    cdef enum VerboseLevel:
        INFO = 1
        WARNING = 2
        ERROR = 3
        NONE = 4

    cdef enum ErrorCode:
        SomethingWentWrong_reportToDeveloper                = 0
        AlignmentNotLoaded                                  = 1
        NoFormatsSpecified                                  = 2
        AlternativeMatrixNotRecognized                      = 3
        ReferenceFileNotLoaded                              = 4
        GapThresholdOutOfRange                              = 5
        GapThresholdNotRecognized                           = 6
        SimilarityThresholdOutOfRange                       = 7
        SimilarityThresholdNotRecognized                    = 8
        ConsistencyThresholdOutOfRange                      = 9
        ConsistencyThresholdNotRecognized                   = 10
        ConservationThresholdOutOfRange                     = 11
        ConservationThresholdNotRecognized                  = 12
        ResidueOverlapOutOfRange                            = 13
        ResidueOverlapNotRecognized                         = 14
        SequencesOverlapOutOfRange                          = 15
        SequencesOverlapNotRecognized                       = 16
        MaxIdentityOutOfRange                               = 17
        MaxIdentityNotRecognized                            = 18
        ClustersValueOutOfRange                             = 19
        ClustersValueNotRecognized                          = 20
        WindowValueOutOfRange                               = 21
        WindowValueNotRecognized                            = 22
        SelectSeqsNotRecognized                             = 23
        SelectColsNotRecognized                             = 24
        GapWindowValueOutOfRange                            = 25
        GapWindowValueNotRecognized                         = 26
        SimilarityWindowValueOutOfRange                     = 27
        SimilarityWindowValueNotRecognized                  = 28
        ConsistencyWindowValueOutOfRange                    = 27
        ConsistencyWindowValueNotRecognized                 = 28
        BlockSizeOutOfRange                                 = 29
        BlockSizeNotRecognized                              = 30
        InFileComparisonStatistics                          = 31
        IncompatibleArguments                               = 32
        SelectSeqsResAndThresholdIncompatibilities          = 33
        SelectSeqsResAndAutomathedMethodsIncompatibilities  = 34
        SelectSeqsResAndWindowIncompatibilities             = 35
        SelectSeqsResAndOverlapIncompatibilites             = 36
        OnlyOneSequencesSelectionMethodAllowed              = 37
        CombinationAmongTrimmingMethods                     = 38
        AutomathicMethodAndBlock                            = 39
        WindowAndArgumentIncompatibilities                  = 40
        CombinationAmongThresholdsMethods                   = 41
        GeneralAndSpecificWindows                           = 42
        StatisticsArgumentIncompatibilities                 = 43
        TrimmingMethodNeeded                                = 44
        ForceFileWithoutCompareDataset                      = 45
        BacktranslationWithoutMainAlignment                 = 46
        NotAligned                                          = 47
        MatrixGivenWithNoMethodToUseIt                      = 48
        SameNameOutput                                      = 49
        SequenceAndResiduesOverlapMutuallyNeeded            = 50
        OutFileNeededWhenPrintingStatistics                 = 51
        AlignmentTypesNotMatching                           = 52
        BlocksizeTooBig                                     = 53
        ParemeterOnlyOnBacktranslation                      = 54
        ProteinAlignmentMustBeAligned                       = 55
        BacktransAlignIsDNA                                 = 56
        ImpossibleToGenerate                                = 57
        ImpossibleToProcessMatrix                           = 58
        SelectOnlyAccepts                                   = 59
        MoreClustersThanSequences                           = 60
        LeftBoundaryBiggerThanRightBoundary                 = 61
        DifferentNumberOfSequencesInCompareset              = 62
        DifferentSeqsNamesInCompareset                      = 63
        CDScontainsProteinSequences                         = 64
        SequenceContainsGap                                 = 65
        SequenceNotMultipleOfThree                          = 66
        SequenceHasStopCodon                                = 67
        SequenceNotPresentInCDS                             = 68
        UnknownCharacter                                    = 69
        SequencesNotSameSize                                = 70
        IncorrectSymbol                                     = 71
        UndefinedSymbol                                     = 72
        ParameterNotFoundOrRepeated                         = 73
        SimilarityMatrixNotCompatibleWindow                 = 74
        PossibleMissmatch                                   = 75
        BracketsMissmatchFound                              = 76
        UnalignedAlignmentToAlignedFormat                   = 77
        CantOpenFile                                        = 78
        FileIsEmpty                                         = 79
        AlignmentFormatNotRecognized                        = 80
        OutputFormatNotRecognized                           = 81
        OnlyOneFormatOnConsoleOutput                        = 82
        AlignmentNotSaved                                   = 83
        VerboseLevelNotRecognized                           = 84
        NeedToSpecifyVerboseLevel                           = 85
        NoReferenceSequenceForContig                        = 86
        SNPoutOfBounds                                      = 87
        NoInputFile                                         = 88
        ComparesetFailedAlignmentMissing                    = 89
        GapWindowTooBig                                     = 90
        SimilarityWindowTooBig                              = 91
        ConsistencyWindowTooBig                             = 92
        WindowTooBig                                        = 93
        AlignmentIsEmpty                                    = 94
        AlignmentTypeIsUnknown                              = 95
        MultipleOutputFormatsSameName                       = 96
        MultipleInputs                                      = 97
        ReferenceNucleotideNotCorresponding                 = 98
        OverwrittingSNP                                     = 99
        MoreDonorsOnLineThanPresented                       = 100
        MinQualityLesserThan0                               = 101
        MinQualityNotRecognized                             = 102
        MinCoverageLesserThan0                              = 103
        MinCoverageNotRecognized                            = 104
        OnlyValidWithVCF                                    = 105
        TriedRenamingOutputPreventOverride                  = 106
        AbsoluteAndRelativeGapThreshold                     = 107
        AbsoluteGapThresholdLessThanZero                    = 108
        AbsoluteGapThresholdBiggerThanNumberOfSequences     = 109
        AbsoluteGapThresholdNotRecognized                   = 110
        MoreThanOneAutomatedMethod                          = 111
        ForceSelectAndInArgumentsProvided                   = 112
        ComparesetAndInArgumentsProvided                    = 113
        NoResidueSequences                                  = 114

    cdef enum WarningCode:
        RemovingOnlyGapsSequence        = 1
        KeepingOnlyGapsSequence         = 2
        SequenceWillBeCut               = 3
        IncludingIndeterminationSymbols = 4
        LessNucleotidesThanExpected     = 5
        HeaderWillBeCut                 = 6
        DonorAlreadyAdded               = 7
        SNPAlreadApplied                = 8
        OverwrittingFile                = 9
        RenamingOutputPreventOverride   = 10

    cdef enum InfoCode:
        CuttingSequence            = 1
        WindowSizeCompareset       = 2
        AddingSNP                  = 3
        RemovingDuplicateSequences = 4

    cdef cppclass reportWrapper:
        bool CanReport

    cdef cppclass reportManager:
        bool IsDebug
        VerboseLevel Level
        map[InfoCode,    const char*] InfoMessages
        map[WarningCode, const char*] WarningMessages
        map[ErrorCode,   const char*] ErrorMessages


cdef extern from "reportsystem.h" nogil:

    cdef reportManager debug
