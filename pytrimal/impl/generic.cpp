#include <limits.h>
#include <stdint.h>
#include <vector>

#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

#include "generic.h"

namespace statistics {
    GenericSimilarity::GenericSimilarity(Alignment* parentAlignment): Similarity(parentAlignment) {}

    bool GenericSimilarity::calculateVectors(bool cutByGap) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool GenericSimilarity::calculateVectors(int *gaps) ");

        // A similarity matrix must be defined. If not, return false
        if (simMatrix == nullptr)
            return false;

        // Calculate the matrix identity in case it's not done before
        if (matrixIdentity == nullptr)
            calculateMatrixIdentity();

        // Create the variable gaps, in case we want to cut by gaps
        int *gaps = nullptr;

        // Retrieve the gaps values in case we want to set to 0 the similarity value
        // in case the gaps value for that column is bigger or equal to 0.8F
        if (cutByGap)
        {
            if (alig->Statistics->gaps == nullptr)
                alig->Statistics->calculateGapStats();
            gaps = alig->Statistics->gaps->getGapsWindow();
        }

        // Initialize the variables used
        int i, j, k;
        float num, den;

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // Q temporal value
        float Q;
        // Temporal chars that will contain the residues to compare by pair.
        int numA, numB;

        // Calculate the maximum number of gaps a column can have to calculate it's
        //      similarity
        float gapThreshold = 0.8F * alig->numberOfResidues;

        // Cache pointers to matrix rows to avoid dereferencing in inner loops
        float* identityRow;
        float* distRow;

        // Create buffers to store column data
        std::vector<char> colnum = std::vector<char>(alig->originalNumberOfSequences);
        std::vector<char> colgap = std::vector<char>(alig->originalNumberOfSequences);

        // For each column calculate the Q value and the MD value using an equation
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            // Set MDK for columns with gaps values bigger or equal to 0.8F
            if (cutByGap && gaps[i] >= gapThreshold) {
                MDK[i] = 0.F;
                continue;
            }

            // Fill the column buffer with the current column and check
            // characters are well-defined with respect to the similarity matrix
            for (j = 0; j < alig->originalNumberOfSequences; j++) {
                char letter = utils::toUpper(alig->sequences[j][i]);
                if ((letter == indet) || (letter == '-')) {
                    colgap[j] = 1;
                } else {
                    colgap[j] = 0;
                    if ((letter < 'A') || (letter > 'Z')) {
                        debug.report(ErrorCode::IncorrectSymbol, new std::string[1]{std::string(1, letter)});
                        return false;
                    } else if (simMatrix->vhash[letter - 'A'] == -1) {
                        debug.report(ErrorCode::UndefinedSymbol, new std::string[1]{std::string(1, letter)});
                        return false;
                    } else {
                        colnum[j] = simMatrix->vhash[letter - 'A'];
                    }
                }
            }

            // For each AAs/Nucleotides' pair in the column we compute its distance
            for (j = 0, num = 0, den = 0; j < alig->originalNumberOfSequences; j++) {
                // We don't compute the distance if the first element is
                // a indeterminate (XN) or a gap (-) element.
                if (colgap[j]) continue;

                // Get the index of the first residue
                // and cache pointers to matrix rows
                numA = colnum[j];
                distRow = simMatrix->distMat[numA];
                identityRow = matrixIdentity[j];

                for (k = j + 1; k < alig->originalNumberOfSequences; k++) {
                    // We don't compute the distance if the second element is
                    //      a indeterminate (XN) or a gap (-) element
                    if (colgap[k]) continue;

                    // Get the index of the second residue and compute 
                    // fraction with identity value for the two pairs and
                    // its distance based on similarity matrix's value.
                    numB = colnum[k];
                    num += identityRow[k] * distRow[numB];
                    den += identityRow[k];
                }
            }

            // If we are processing a column with only one AA/nucleotide, MDK = 0
            if (den == 0)
                MDK[i] = 0;
            else
            {
                Q = num / den;
                // If the MDK value is more than 1, we normalized this value to 1.
                //      Only numbers higher than 0 yield exponents higher than 1
                //      Using this we can test if the result is going to be higher than 1.
                //      And thus, prevent calculating the exp.
                // Take in mind that the Q is negative, so we must test if Q is LESSER
                //      than one, not bigger.
                if (Q < 0)
                    MDK[i] = 1.F;
                else
                    MDK[i] = exp(-Q);
            }
        }

        for (i = 0; i < alig->originalNumberOfSequences; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
        matrixIdentity = nullptr;

        return true;
    }
}

GenericCleaner::GenericCleaner(Alignment* parent): Cleaner(parent) {
    hits = new uint32_t[alig->originalNumberOfResidues];
}

GenericCleaner::~GenericCleaner() {
    delete[] hits;
}

bool GenericCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool GenericCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) ");

    // abort if there is not output vector to write to
    if (spuriousVector == nullptr)
        return false;

    // compute number of sequences from overlap threshold
    uint32_t ovrlap  = uint32_t(ceil(overlap * float(alig->originalNumberOfSequences - 1)));

    // Depending on alignment type, indetermination symbol will be one or other
    char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

    // for each sequence in the alignment, computes its overlap
    for (int i = 0; i < alig->originalNumberOfSequences; i++) {

        // reset hits count
        memset(&hits[0], 0, alig->originalNumberOfResidues*sizeof(uint32_t));

        // compare sequence to other sequences for every position
        for (int j = 0; j < alig->originalNumberOfSequences; j++) {

            // don't compare sequence to itself
            if (j == i)
                continue;

            const char* datai = alig->sequences[i].data();
            const char* dataj = alig->sequences[j].data();

            // process the tail elements when there remain less than
            // can be fitted in a SIMD vector
            for (int k = 0; k < alig->originalNumberOfResidues; k++) {
                int nongapi = (datai[k] != indet) && (datai[k] != '-');
                int nongapj = (dataj[k] != indet) && (dataj[k] != '-');
                hits[k] += ((nongapi && nongapj) || (datai[k] == dataj[k]));
            }
        }

        // compute number of good positions in for sequence i
        uint32_t seqValue = 0;
        for (int k = 0; k < alig->originalNumberOfResidues; k++)
            if (hits[k] >= ovrlap)
                seqValue++;

        // compute overlap of current sequence as the fraction of columns
        // above overlap threshold
        spuriousVector[i] = ((float) seqValue / alig->originalNumberOfResidues);
    }

    // If there is not problem in the method, return true
    return true;
}
