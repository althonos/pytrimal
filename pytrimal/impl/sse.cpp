#include <immintrin.h>
#include <limits.h>
#include <stdint.h>

#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

#include "sse.h"

static inline uint32_t _mm_hsum_epi8(__m128i a) {
    __m128i vsum = _mm_sad_epu8(a, _mm_setzero_si128());
    return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
}

namespace statistics {
    SSESimilarity::SSESimilarity(Alignment* parentAlignment): Similarity(parentAlignment) {
        ascii_vhash = std::vector<char>(UCHAR_MAX, -1);
        column = std::string(parentAlignment->originalNumberOfSequences, 0);
        colgap = std::vector<char>(parentAlignment->originalNumberOfSequences);
    }

    void SSESimilarity::calculateMatrixIdentity() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void SSESimilarity::calculateMatrixIdentity() ");

        // We don't want to calculate the matrix identity
        // if it has been previously calculated
        if (matrixIdentity != nullptr)
            return;

        // declare indices
        int i, j, k, l;

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

                // run with unrolled loops of UCHAR_MAX iterations first
                for (k = 0; k + ((int) sizeof(__m128i))*UCHAR_MAX < alig->originalNumberOfResidues;) {
                    // unroll the internal loop
                    for (l = 0; l < UCHAR_MAX; l++, k += sizeof(__m128i)) {
                        // load data for the sequences
                        seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                        seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
                        // find which sequence characters are gap or indet
                        gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                        gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                        // find which sequence characters are equal
                        eq    = _mm_cmpeq_epi8(seqi, seqj);
                        // update counters
                        sum_acc = _mm_add_epi8(sum_acc, _mm_and_si128(eq, _mm_andnot_si128( _mm_or_si128(gapsi, gapsj), ONES)));
                        len_acc = _mm_add_epi8(len_acc,                   _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
                    }
                    // merge accumulators
                    sum    += _mm_hsum_epi8(sum_acc);
                    length += _mm_hsum_epi8(len_acc);
                    sum_acc = _mm_setzero_si128();
                    len_acc = _mm_setzero_si128();
                }

                // run remaining iterations in SIMD while possible
                for (; k + ((int) sizeof(__m128i)) < alig->originalNumberOfResidues; k += sizeof(__m128i)) {
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
                }

                // merge accumulators
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

    bool SSESimilarity::calculateVectors(bool cutByGap) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool SSESimilarity::calculateVectors(int *gaps) ");

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
        char chA, chB;
        int numA, numB;

        // Calculate the maximum number of gaps a column can have to calculate it's
        //      similarity
        float gapThreshold = 0.8F * alig->numberOfResidues;

        // For each column calculate the Q value and the MD value using an equation
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            // Set MDK for columns with gaps values bigger or equal to 0.8F
            if (cutByGap && (gaps[i] >= gapThreshold)) {
                MDK[i] = 0.F;
                continue;
            }

            // Fill the column data with the current column and check characters
            // are well-defined with respect to the similarity matrix
            for (j = 0; j < alig->originalNumberOfSequences; j++) {
                column[j] = chA = utils::toUpper(alig->sequences[j][i]);
                if ((chA == indet) || (chA == '-')) {
                    colgap[j] = 1;
                } else {
                    colgap[j] = 0;
                    if ((chA < 'A') || (chA > 'Z')) {
                        debug.report(ErrorCode::IncorrectSymbol, new std::string[1]{std::string(1, chA)});
                        return false;
                    } else if (ascii_vhash[chA] == -1) {
                        debug.report(ErrorCode::UndefinedSymbol, new std::string[1]{std::string(1, chA)});
                        return false;
                    }
                }
            }

            // For each AAs/Nucleotides' pair in the column we compute its distance
            for (j = 0, num = 0, den = 0; j < alig->originalNumberOfSequences; j++) {
                // We don't compute the distance if the first element is
                // a indeterminate (XN) or a gap (-) element.
                if (colgap[j]) continue;

                // Calculate the upper value of the residue,
                //      to use in simMatrix->getDistance
                // This is faster than calculating the upper on that method
                //      as this is done before entering the loop
                // Doing this before checking if the element is indeterminate or gap
                //      allows to check if the indetermination is not capitalized
                chA = column[j];
                // Search the first character position
                numA = ascii_vhash[chA];

                for (k = j + 1; k < alig->originalNumberOfSequences; k++) {
                    // We don't compute the distance if the second element is
                    //      a indeterminate (XN) or a gap (-) element
                    if (colgap[k]) continue;

                    // We calculate the upper value of the residue,
                    //      to use in simMatrix->getDistance
                    // This is equally faster as if it was done inside the method
                    //      but to prevent errors, the method doesn't 'upper'
                    //      the given chars.
                    // Doing this before checking if the element is indeterminate or gap
                    //      allows to check if the indetermination is not capitalized
                    chB = column[k];
                    // Search the second character position
                    numB = ascii_vhash[chB];

                    // We use the identity value for the two pairs and
                    //      its distance based on similarity matrix's value.
                    num += matrixIdentity[j][k] * simMatrix->distMat[numA][numB];
                    den += matrixIdentity[j][k];
                }
            }

            // If we are processing a column with only one AA/nucleotide, MDK = 0
            if (den == 0) {
                MDK[i] = 0;
            } else {
                Q = num / den;
                // If the MDK value is more than 1, we normalized this value to 1.
                //      Only numbers higher than 0 yield exponents higher than 1
                //      Using this we can test if the result is going to be higher than 1.
                //      And thus, prevent calculating the exp.
                // Take in mind that the Q is negative, so we must test if Q is LESSER
                //      than one, not bigger.
                MDK[i] = (Q < 0) ? 1.F : exp(-Q);;
            }
        }

        for (i = 0; i < alig->originalNumberOfSequences; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
        matrixIdentity = nullptr;

        return true;
    }

    bool SSESimilarity::setSimilarityMatrix(similarityMatrix *sm) {
        if (sm != nullptr) {
            for (char x = 'A'; x <= 'Z'; x++)
                ascii_vhash[x] = sm->vhash[x - 'A'];
        }
        return Similarity::setSimilarityMatrix(sm);
    }
}

SSECleaner::SSECleaner(Alignment* parent): Cleaner(parent) {
    // create an index for residues to skip
    skipResidues = (unsigned char*) new unsigned char[alig->originalNumberOfResidues];
    for(int i = 0; i < alig->originalNumberOfResidues; i++) {
        skipResidues[i] = alig->saveResidues[i] == -1 ? 0xFF : 0;
    }
}

SSECleaner::~SSECleaner() {
    delete[] skipResidues;
}

void SSECleaner::calculateSeqIdentity() {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void SSECleaner::calculateSeqIdentity(void) ");

  // declare indices
  int i, j, k, l;

  // Depending on alignment type, indetermination symbol will be one or other
  char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

  // prepare constant SIMD vectors
  const __m128i allindet = _mm_set1_epi8(indet);
  const __m128i allgap   = _mm_set1_epi8('-');
  const __m128i ONES     = _mm_set1_epi8(1);

  // Create identities matrix to store identities scores
  alig->identities = new float*[alig->originalNumberOfSequences];
  for(i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;
      alig->identities[i] = new float[alig->originalNumberOfSequences];
      alig->identities[i][i] = 0;
  }

  // For each seq, compute its identity score against the others in the MSA
  for (i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;

      // Compute identity scores for the current sequence against the rest
      for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
          if (alig->saveSequences[j] == -1) continue;

          const char* datai = alig->sequences[i].data();
          const char* dataj = alig->sequences[j].data();

          __m128i seqi;
          __m128i seqj;
          __m128i skip;
          __m128i eq;
          __m128i gapsi;
          __m128i gapsj;
          __m128i mask;
          __m128i dst_acc = _mm_setzero_si128();
          __m128i hit_acc = _mm_setzero_si128();

          int hit = 0;
          int dst = 0;

          // run with unrolled loops of UCHAR_MAX iterations first
          for (k = 0; k + ((int) sizeof(__m128i))*UCHAR_MAX < alig->originalNumberOfResidues;) {
              for (l = 0; l < UCHAR_MAX; l++, k += sizeof(__m128i)) {
                  // load data for the sequences
                  seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                  seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
                  skip = _mm_loadu_si128( (const __m128i*) (&skipResidues[k]));
                  eq = _mm_cmpeq_epi8(seqi, seqj);
                  // find which sequence characters are gap or indet
                  gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                  gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                  // find position where not both characters are gap
                  mask = _mm_andnot_si128(skip, _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
                  // update counters
                  dst_acc = _mm_add_epi8(dst_acc,                   mask);
                  hit_acc = _mm_add_epi8(hit_acc, _mm_and_si128(eq, mask));
              }
              // merge accumulators
              dst += _mm_hsum_epi8(dst_acc);
              hit += _mm_hsum_epi8(hit_acc);
              dst_acc = _mm_setzero_si128();
              hit_acc = _mm_setzero_si128();
          }

          // run remaining iterations in SIMD while possible
          for (; k + ((int) sizeof(__m128i)) < alig->originalNumberOfResidues; k += sizeof(__m128i)) {
              // load data for the sequences; load is unaligned because
              // string data is not guaranteed to be aligned.
              seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
              seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
              skip = _mm_loadu_si128( (const __m128i*) (&skipResidues[k]));
              eq = _mm_cmpeq_epi8(seqi, seqj);
              // find which sequence characters are either a gap or
              // an indeterminate character
              gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
              gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
              // find position where not both characters are gap
              mask = _mm_andnot_si128(skip, _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
              // update counters: update dst if any of the two sequence
              // characters is non-gap/indeterminate; update length if
              // any of the characters is non-gappy/indeterminate.
              dst_acc = _mm_add_epi8(dst_acc,                   mask);
              hit_acc = _mm_add_epi8(hit_acc, _mm_and_si128(eq, mask));
          }

          // update counters after last loop
          hit += _mm_hsum_epi8(hit_acc);
          dst += _mm_hsum_epi8(dst_acc);

          // process the tail elements when there remain less than
          // can be fitted in a SIMD vector
          for (; k < alig->originalNumberOfResidues; k++) {
              int gapi = (datai[k] == indet) || (datai[k] == '-');
              int gapj = (dataj[k] == indet) || (dataj[k] == '-');
              dst += (!(gapi && gapj)) && (!skipResidues[k]);
              hit += (!(gapi && gapj)) && (!skipResidues[k]) && (datai[k] == dataj[k]);
          }

          if (dst == 0) {
              debug.report(ErrorCode::NoResidueSequences,
                  new std::string[2]
                      {
                          alig->seqsName[i],
                          alig->seqsName[j]
                      });
              alig->identities[i][j] = 0;
          } else {
              // Identity score between two sequences is the ratio of identical residues
              // by the total length (common and no-common residues) among them
              alig->identities[i][j] = (float) hit / dst;
          }

          alig->identities[j][i] = alig->identities[i][j];
      }
  }
}
