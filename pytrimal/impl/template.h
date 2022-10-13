#include <cstdlib>

#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "Statistics/similarityMatrix.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

namespace simd {

    template<class T, class Vector>
    T* aligned_array(size_t n) {
        const size_t mask = Vector::SIZE - 1;
        const size_t size = (n * sizeof(T) + mask) & (~mask);
        return static_cast<T*>(aligned_alloc(Vector::SIZE, size));
    }

    template<class Vector>
    inline void calculateMatrixIdentity(const Alignment* alig, float** matrixIdentity) {
        // declare indices
        int i, j, k, l;

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // prepare constant SIMD vectors
        const Vector allindet = Vector(indet);
        const Vector ALLGAP   = Vector('-');
        const Vector ONES     = Vector(1);

        // For each sequences' pair, compare identity
        for (i = 0; i < alig->originalNumberOfSequences; i++) {

            const uint8_t* datai = reinterpret_cast<const uint8_t*>(alig->sequences[i].data());

            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {

                const uint8_t* dataj = reinterpret_cast<const uint8_t*>(alig->sequences[j].data());

                Vector len_acc = Vector();
                Vector sum_acc = Vector();

                uint32_t sum    = 0;
                uint32_t length = 0;

                // run with unrolled loops of UCHAR_MAX iterations first
                for (k = 0; ((int) (k + Vector::LANES*UCHAR_MAX)) < alig->originalNumberOfResidues;) {
                    // unroll the internal loop
                    for (l = 0; l < UCHAR_MAX; l++, k += Vector::LANES) {
                        // load data for the sequences
                        Vector seqi = Vector::load(&datai[k]);
                        Vector seqj = Vector::load(&dataj[k]);
                        // find which sequence characters are gap or indet
                        Vector gapsi = (seqi == ALLGAP) | (seqi == allindet);
                        Vector gapsj = (seqj == ALLGAP) | (seqj == allindet);
                        // find which sequence characters are equal
                        Vector eq    = (seqi == seqj);
                        // update counters
                        sum_acc += eq & ONES.andnot(gapsi | gapsj);
                        len_acc += ONES.andnot(gapsi & gapsj);
                    }
                    // merge accumulators
                    sum    += sum_acc.sum();
                    length += len_acc.sum();
                    sum_acc.clear();
                    len_acc.clear();
                }

                // run remaining iterations in SIMD while possible
                for (; ((int) (k + Vector::LANES)) < alig->originalNumberOfResidues; k += Vector::LANES) {
                    // load data for the sequences
                    Vector seqi = Vector::load(&datai[k]);
                    Vector seqj = Vector::load(&dataj[k]);
                    // find which sequence characters are gap or indet
                    Vector gapsi = (seqi == ALLGAP) | (seqi == allindet);
                    Vector gapsj = (seqj == ALLGAP) | (seqj == allindet);
                    // find which sequence characters are equal
                    Vector eq    = (seqi == seqj);
                    // update counters
                    sum_acc += eq & ONES.andnot(gapsi | gapsj);
                    len_acc += ONES.andnot(gapsi & gapsj);
                }

                // // merge accumulators
                sum    += sum_acc.sum();
                length += len_acc.sum();

                // process the tail elements when there remain less than
                // can be fitted in a SIMD vector
                for (; k < alig->originalNumberOfResidues; k++) {
                    int gapi = (datai[k] == '-') || (datai[k] == indet);
                    int gapj = (dataj[k] == '-') || (dataj[k] == indet);
                    sum    += (!gapi) && (!gapj) && (datai[k] == dataj[k]);
                    length += (!gapi) || (!gapj);
                }

                // Calculate the value of matrix idn for columns j and i
                matrixIdentity[i][j] = matrixIdentity[j][i] = (1.0F - ((float)sum / length));
            }
        }
    }

    template<class Vector>
    inline bool calculateSpuriousVector(const Alignment* alig, const float overlap, float *spuriousVector) {
        // compute number of sequences from overlap threshold
        uint32_t  ovrlap  = uint32_t(ceil(overlap * float(alig->originalNumberOfSequences - 1)));

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

        // prepare constant SIMD vectors
        const Vector allindet = Vector(indet);
        const Vector ALLGAP   = Vector('-');
        const Vector ONES     = Vector(1);

        // allocate aligned memory for faster SIMD loads
        uint32_t* hits    = aligned_array<uint32_t, Vector>(alig->originalNumberOfResidues);
        uint8_t*  hits_u8 = aligned_array<uint8_t, Vector>(alig->originalNumberOfResidues);

        // for each sequence in the alignment, computes its overlap
        for (int i = 0; i < alig->originalNumberOfSequences; i++) {

            // reset hits count
            memset(&hits[0],    0, alig->originalNumberOfResidues*sizeof(uint32_t));
            memset(&hits_u8[0], 0, alig->originalNumberOfResidues*sizeof(uint8_t));

            const uint8_t* datai = reinterpret_cast<const uint8_t*>(alig->sequences[i].data());

            // compare sequence to other sequences for every position
            for (int j = 0; j < alig->originalNumberOfSequences; j++) {

                // don't compare sequence to itself
                if (j == i)
                    continue;

                const uint8_t* dataj = reinterpret_cast<const uint8_t*>(alig->sequences[j].data());

                int k = 0;

                // run iterations in SIMD while possible
                for (; ((int) (k + Vector::LANES)) <= alig->originalNumberOfResidues; k += Vector::LANES) {
                    // load data for the sequences
                    const Vector seqi = Vector::loadu(&datai[k]);
                    const Vector seqj = Vector::loadu(&dataj[k]);
                    // find which sequence characters are gap or indet
                    const Vector gapsi = (seqi == ALLGAP) | (seqi == allindet);
                    const Vector gapsj = (seqj == ALLGAP) | (seqj == allindet);
                    const Vector gaps  = !(gapsi | gapsj);
                    // find which sequence characters match
                    const Vector eq = (seqi == seqj);
                    // find position where either not both characters are gap, or they are equal
                    const Vector n = (eq | gaps) & ONES;
                    // update counters
                    Vector hit = Vector::load(&hits_u8[k]);
                    hit += n;
                    hit.store(&hits_u8[k]);
                }

                // process the tail elements when there remain less than
                // can be fitted in a SIMD vector
                for (; k < alig->originalNumberOfResidues; k++) {
                    int nongapi = (datai[k] != indet) && (datai[k] != '-');
                    int nongapj = (dataj[k] != indet) && (dataj[k] != '-');
                    hits_u8[k] += ((nongapi && nongapj) || (datai[k] == dataj[k]));
                }

                // we can process up to UCHAR_MAX sequences, otherwise hits_u8[k]
                // may overflow, so every UCHAR_MAX iterations we transfer the
                // partial hit counts from `hits_u8` to `hits`
                if ((j % UCHAR_MAX) == 0) {
                    for (k = 0; k < alig->originalNumberOfResidues; k++) hits[k] += hits_u8[k];
                    memset(hits_u8, 0, alig->originalNumberOfResidues*sizeof(uint8_t));
                }
            }

            // update counters after last loop
            for (int k = 0; k < alig->originalNumberOfResidues; k++) hits[k] += hits_u8[k];

            // compute number of good positions in for sequence i
            uint32_t seqValue = 0;
            for (int k = 0; k < alig->originalNumberOfResidues; k++)
                if (hits[k] >= ovrlap)
                    seqValue++;

            // compute overlap of current sequence as the fraction of columns
            // above overlap threshold
            spuriousVector[i] = ((float) seqValue / alig->originalNumberOfResidues);
        }

        // free allocated memory
        free(hits);
        free(hits_u8);

        // If there is not problem in the method, return true
        return true;
    }

    template<class Vector>
    inline void calculateSeqIdentity(const Alignment* alig, float** identities) {
            
        // declare indices
        int i, j, k, l;

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

        // prepare constant SIMD vectors
        const Vector allindet = Vector(indet);
        const Vector ALLGAP   = Vector('-');
        const Vector ONES     = Vector(1);

        // create an index of residues to skip
        uint8_t* skipResidues = aligned_array<uint8_t, Vector>(alig->originalNumberOfResidues); //ALIGNED_ALLOC(alig->originalNumberOfResidues, uint8_t);
        for(k = 0; k < alig->originalNumberOfResidues; k++) {
            skipResidues[k] = alig->saveResidues[k] == -1 ? 0xFF : 0;
        }

        // For each seq, compute its identity score against the others in the MSA
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            if (alig->saveSequences[i] == -1) continue;

            const uint8_t* datai = reinterpret_cast<const uint8_t*>(alig->sequences[i].data());


            // Compute identity scores for the current sequence against the rest
            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
                if (alig->saveSequences[j] == -1) continue;

                const uint8_t* dataj = reinterpret_cast<const uint8_t*>(alig->sequences[j].data());

                Vector dst_acc = Vector();
                Vector hit_acc = Vector();

                int hit = 0;
                int dst = 0;

                // run with unrolled loops of UCHAR_MAX iterations first
                for (k = 0; ((int) (k + Vector::LANES*UCHAR_MAX)) < alig->originalNumberOfResidues;) {
                    for (l = 0; l < UCHAR_MAX; l++, k += Vector::LANES) {
                        // load data for the sequences
                        Vector seqi = Vector::loadu( (&datai[k]));
                        Vector seqj = Vector::loadu( (&dataj[k]));
                        Vector skip = Vector::load(  (&skipResidues[k]));
                        Vector eq = (seqi == seqj);
                        // find which sequence characters are gap or indet
                        Vector gapsi = ((seqi == ALLGAP) | (seqi == allindet));
                        Vector gapsj = ((seqj == ALLGAP) | (seqj == allindet));
                        // find position where not both characters are gap
                        Vector mask = ONES.andnot(gapsi & gapsj).andnot(skip);
                        // update counters
                        dst_acc += mask;
                        hit_acc += (eq & mask);
                    }
                    // merge accumulators
                    dst += dst_acc.sum();
                    hit += hit_acc.sum();
                    dst_acc.clear();
                    hit_acc.clear();
                }

                // run remaining iterations in SIMD while possible
                for (; ((int) (k + Vector::LANES)) < alig->originalNumberOfResidues; k += Vector::LANES) {
                    // load data for the sequences; load is unaligned because
                    // string data is not guaranteed to be aligned.
                    Vector seqi = Vector::loadu( (&datai[k]));
                    Vector seqj = Vector::loadu( (&dataj[k]));
                    Vector skip = Vector::load(  (&skipResidues[k]));
                    Vector eq = (seqi == seqj);
                    // find which sequence characters are gap or indet
                        Vector gapsi = ((seqi == ALLGAP) | (seqi == allindet));
                        Vector gapsj = ((seqj == ALLGAP) | (seqj == allindet));
                        // find position where not both characters are gap
                        Vector mask = ONES.andnot(gapsi & gapsj).andnot(skip);
                        // update counters
                        dst_acc += mask;
                        hit_acc += (eq & mask);
                }

                // update counters after last loop
                hit += hit_acc.sum();
                dst += dst_acc.sum();

                // process the tail elements when there remain less than
                // can be fitted in a SIMD vector
                for (; k < alig->originalNumberOfResidues; k++) {
                    int gapi = (datai[k] == indet) || (datai[k] == '-');
                    int gapj = (dataj[k] == indet) || (dataj[k] == '-');
                    dst += (!(gapi && gapj)) && (!skipResidues[k]);
                    hit += (!(gapi && gapj)) && (!skipResidues[k]) && (datai[k] == dataj[k]);
                }

                if (dst == 0) {
                    debug.report(
                            ErrorCode::NoResidueSequences,
                        new std::string[2] { alig->seqsName[i], alig->seqsName[j] }
                        );
                    identities[i][j] = 0;
                } else {
                    // Identity score between two sequences is the ratio of identical residues
                    // by the total length (common and no-common residues) among them
                    alig->identities[i][j] = (float) hit / dst;
                }

                identities[j][i] = identities[i][j];
            }
        }

        // free allocated memory
        free(skipResidues);

    }

    template<class Vector>
    inline void calculateGapVectors(const Alignment* alig, int* gapsInColumn) {
        int i, j;
        const __m128i ALLGAP = _mm_set1_epi8('-');
        const __m128i ONES   = _mm_set1_epi8(1);

        // use temporary buffer for storing 8-bit partial sums
        uint8_t* gapsInColumn_u8 = aligned_array<uint8_t, Vector>(alig->originalNumberOfResidues);
        memset(gapsInColumn,    0, sizeof(int)     * alig->originalNumberOfResidues);
        memset(gapsInColumn_u8, 0, sizeof(uint8_t) * alig->originalNumberOfResidues);

        // count gaps per column
        for (j = 0; j < alig->originalNumberOfSequences; j++) {
            // skip sequences not retained in alignment
            if (alig->saveSequences[j] == -1)
                continue;
            // process the whole sequence, 16 lanes at a time
            const char* data = alig->sequences[j].data();
            for (i = 0; ((int) (i + Vector::LANES)) < alig->originalNumberOfResidues; i += Vector::LANES) {
                __m128i letters = _mm_loadu_si128((const __m128i*) &data[i]);
                __m128i counts  = _mm_load_si128((const __m128i*) &gapsInColumn_u8[i]);
                __m128i gaps    = _mm_and_si128(ONES, _mm_cmpeq_epi8(letters, ALLGAP));
                __m128i updated = _mm_add_epi8(counts, gaps);
                _mm_store_si128((__m128i*) &gapsInColumn_u8[i], updated);
            }
            // count the remaining gap elements without SIMD
            for (; i < alig->originalNumberOfResidues; i++)
                if (data[i] == '-') gapsInColumn_u8[i]++;
            // every UCHAR_MAX iterations the accumulator may overflow, so the
            // temporary counts are moved into the final counter array, and the
            // accumulator is reset
            if (j % UCHAR_MAX == 0) {
                for (i = 0; i < alig->originalNumberOfResidues; i++)
                    gapsInColumn[i] += gapsInColumn_u8[i];
                memset(gapsInColumn_u8, 0, sizeof(uint8_t) * alig->originalNumberOfResidues);
            }
        }
        // collect the remaining partial counts into the final counter array
        for (i = 0; i < alig->originalNumberOfResidues; i++)
            gapsInColumn[i] += gapsInColumn_u8[i];

        // free temporary buffer
        free(gapsInColumn_u8);
    }

    template<class Vector>
    inline bool calculateSimilarityVectors(const Alignment* alig, const float** matrixIdentity, const statistics::similarityMatrix* simMatrix, const int* gaps, float* MDK) {
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
        const float* identityRow;
        const float* distRow;

        // Create buffers to store column data
        std::vector<char> colnum = std::vector<char>(alig->originalNumberOfSequences);
        std::vector<char> colgap = std::vector<char>(alig->originalNumberOfSequences);

        // For each column calculate the Q value and the MD value using an equation
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            // Set MDK for columns with gaps values bigger or equal to threshold
            if ((gaps != nullptr) && gaps[i] >= gapThreshold) {
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

        return true;
    }
}