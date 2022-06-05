#include <immintrin.h>
#include <stdint.h>

#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"

#include "sse.h"

void print_cvec(__m128i v) {
    union { __m128i vec; unsigned char text[16]; } x;
    x.vec = v;

    for (int i = 0; i < 16; i++)
        printf("%c ", x.text[i]);
    printf("\n");
}

void print_vec(__m128i v) {
    union { __m128i vec; unsigned char text[16]; } x;
    x.vec = v;

    for (int i = 0; i < 16; i++)
        printf("%3u ", x.text[i]);
    printf("\n");
}

static inline uint64_t _mm_hsum_epi8(__m128i a) {
  const __m128i ZEROS = _mm_setzero_si128();

  __m128i x = _mm_unpackhi_epi8(a, ZEROS); // [ 0 a0  0 a1  0  a2  0  a3  0  a4  0  a5  0  a6  0  a7]
  __m128i y = _mm_unpacklo_epi8(a, ZEROS); // [ 0 a8  0 a9  0 a10  0 a11  0 a12  0 a13  0 a14  0 a15]
  __m128i z = _mm_add_epi16(x, y);         // [a0+a8 a1+a9 a2+a10 a3+a11 a4+a12 a5+a13 a6+a14 a7+a15]
                                           // [  X0     X1     X2     X3     X4     X5     X6     X7]

  __m128i h1 = _mm_add_epi16(z,  _mm_shuffle_epi32(z,  0b10110001));  // [ X0+X2 X1+X3 X2+X0 X3+X1 X4+X6 X5+X7 X6+X4 X7+X5 ]
  __m128i h2 = _mm_add_epi16(h1, _mm_shuffle_epi32(h1, 0b00011011));  // [ X0+X2+X6+X4 X1+X3+X5+X7 ....]
  return _mm_extract_epi16(h2, 0) + _mm_extract_epi16(h2, 1);         // (X0+X2+X6+X4) + (X1+X3+X5+X7)
}

namespace statistics {

    SSESimilarity::SSESimilarity(Alignment* parentAlignment): Similarity(parentAlignment) {}

    void SSESimilarity::calculateMatrixIdentity() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void SSESimilarity::calculateMatrixIdentity() ");

        // We don't want to calculate the matrix identity
        // if it has been previously calculated
        if (matrixIdentity != nullptr)
            return;

        // Allocate temporal variables
        int i, j, k;

        // Allocate memory for the matrix identity
        matrixIdentity = new float *[alig->originalNumberOfSequences];
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            matrixIdentity[i] = new float[alig->originalNumberOfSequences];
        }

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // prepare constant SIMD vectors
        const __m128i allindet = _mm_set1_epi8(indet);
        const __m128i allgap   = _mm_set1_epi8('-');
        const __m128i ONES     = _mm_set1_epi8(1);

        // For each sequences' pair, compare identity
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {

                const char* datai = alig->sequences[i].data();
                const char* dataj = alig->sequences[j].data();

                __m128i seqi;
                __m128i seqj;
                __m128i eq;
                __m128i gapsi;
                __m128i gapsj;
                __m128i len_acc = _mm_setzero_si128();
                __m128i sum_acc = _mm_setzero_si128();

                int sum = 0;
                int length = 0;

                for (k = 0; k + sizeof(__m128i) < alig->originalNumberOfResidues; k += sizeof(__m128i)) {
                    // load data for the sequences; load is unaligned because
                    // string data is not guaranteed to be aligned.
                    seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                    seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));

                    // find which sequence characters are either a gap or
                    // an indeterminate character
                    gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                    gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                    // find which sequence characters are equal
                    eq    = _mm_cmpeq_epi8(seqi, seqj);

                    // update counters: update sum if both sequence characters
                    // are non-gap/indeterminate and equal; update length if
                    // any of the characters is non-gappy/indeterminate.
                    sum_acc = _mm_add_epi8(sum_acc, _mm_and_si128(eq, _mm_andnot_si128( _mm_or_si128(gapsi, gapsj), ONES)));
                    len_acc = _mm_add_epi8(len_acc,                   _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));

                    // at least 255 iterations can be done until the accumulators
                    // overflow; when we reach iteration 255, we move data from
                    // the lane accumulators to the main counter variables
                    if (k % (255 * sizeof(__m128i))) {
                        sum    += _mm_hsum_epi8(sum_acc);
                        length += _mm_hsum_epi8(len_acc);
                        sum_acc = _mm_setzero_si128();
                        len_acc = _mm_setzero_si128();
                    }
                }

                // update counters after last loop
                sum    += _mm_hsum_epi8(sum_acc);
                length += _mm_hsum_epi8(len_acc);

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
}
