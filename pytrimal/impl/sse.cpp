#include <immintrin.h>
#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

#include "sse.h"
#include "template.h"

class SSEVector {
private:
    __m128i vector;
    inline SSEVector(__m128i vec): vector(vec) {}
public:
    const static size_t LANES = 16;
    const static size_t SIZE = sizeof(__m128i);

    inline SSEVector(): vector(_mm_setzero_si128()) {}
    inline SSEVector(const int8_t value): vector(_mm_set1_epi8(value)) {} 

    inline static SSEVector load(const uint8_t* data) {
        return SSEVector(_mm_load_si128((const __m128i*) data));
    }

    inline static SSEVector loadu(const uint8_t* data) {
        return SSEVector(_mm_loadu_si128((const __m128i*) data));
    }

    inline void store(uint8_t* data) const {
        _mm_store_si128((__m128i*) data, vector);
    }

    inline void storeu(uint8_t* data) const {
        _mm_storeu_si128((__m128i*) data, vector);
    }

    inline SSEVector& operator+=(const SSEVector& rhs) {
        vector = _mm_add_epi8(vector, rhs.vector);
        return *this;
    }

    inline SSEVector operator==(const SSEVector& rhs) const {
        return SSEVector(_mm_cmpeq_epi8(vector, rhs.vector));
    }

    inline SSEVector operator&(const SSEVector& rhs) const {
        return SSEVector(_mm_and_si128(vector, rhs.vector));
    }

    inline SSEVector operator|(const SSEVector& rhs) const {
        return SSEVector(_mm_or_si128(vector, rhs.vector));
    }

    inline SSEVector operator!() const {
        return SSEVector(_mm_andnot_si128(vector, _mm_set1_epi8(0xFF)));
    }

    inline SSEVector andnot(const SSEVector& rhs) const {
        return SSEVector(_mm_andnot_si128(rhs.vector, vector));
    }

    inline uint16_t sum() const {
        __m128i vsum = _mm_sad_epu8(vector, _mm_setzero_si128());
        return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
    }

    inline void clear() {
        vector = _mm_setzero_si128();
    }

};

namespace statistics {
    void SSESimilarity::calculateMatrixIdentity() {
        // create a timer for this function
        StartTiming("void SSESimilarity::calculateMatrixIdentity() ");
        // abort if identity matrix computation was already done
        if (matrixIdentity != nullptr)
            return;
        // Allocate memory for the matrix identity
        matrixIdentity = new float *[alig->originalNumberOfSequences];
        for (int i = 0; i < alig->originalNumberOfSequences; i++) {
            matrixIdentity[i] = new float[alig->originalNumberOfSequences];
        }
        // Run SIMD code with SSE
        simd::calculateMatrixIdentity<SSEVector>(alig, matrixIdentity);
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

        if (!simd::calculateSimilarityVectors<SSEVector>(alig, const_cast<const float**>(matrixIdentity), simMatrix, gaps, MDK)) 
            return false;

        for (int i = 0; i < alig->originalNumberOfSequences; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
        matrixIdentity = nullptr;

        return true;
    }

    void SSEGaps::CalculateVectors() {
        // create a timer for this function
        StartTiming("bool SSEGaps::CalculateVectors(int *gaps) ");
        // calculate gaps in SIMD with SSE
        simd::calculateGapVectors<SSEVector>(alig, gapsInColumn);
        // build histogram and find largest number of gaps
        for (int i = 0; i < alig->originalNumberOfResidues; i++) {
            totalGaps += gapsInColumn[i];
            numColumnsWithGaps[gapsInColumn[i]]++;
            if (gapsInColumn[i] > maxGaps)
                maxGaps = gapsInColumn[i];
        }
    }
}

void SSECleaner::calculateSeqIdentity() {
  // create a timer for this function
  StartTiming("void SSECleaner::calculateSeqIdentity(void) ");
  // create identities matrix to store identities scores
  alig->identities = new float*[alig->originalNumberOfSequences];
  for(int i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;
      alig->identities[i] = new float[alig->originalNumberOfSequences];
      alig->identities[i][i] = 0;
  }
  // Run SIMD code with SSE
  simd::calculateSeqIdentity<SSEVector>(alig, alig->identities);
}

bool SSECleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
    // create a timer for this function
    StartTiming("bool SSECleaner::calculateSpuriousVector(float overlap, float *spuriousVector) ");
    // abort if there is not output vector to write to
    if (spuriousVector == nullptr)
        return false;
    // Run SIMD code with SSE
    return simd::calculateSpuriousVector<SSEVector>(alig, overlap, spuriousVector);
}
