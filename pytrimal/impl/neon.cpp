#include <arm_neon.h>
#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

#include "neon.h"

class NEONVector {
private:
    uint8x16_t vector;
    inline NEONVector(uint8x16_t vec): vector(vec) {}
public:
    const static size_t LANES = 16;
    const static size_t SIZE = sizeof(uint8x16_t);

    inline NEONVector(): vector(vdupq_n_u8(0)) {}
    inline NEONVector(const int8_t value): vector(vdupq_n_u8(value)) {} 

    inline static NEONVector load(const uint8_t* data) {
        return NEONVector(vld1q_u8(data));
    }

    inline static NEONVector loadu(const uint8_t* data) {
        return NEONVector(vld1q_u8(data));
    }

    inline void store(uint8_t* data) const {
        vst1q_u8(data, vector);
    }

    inline void storeu(uint8_t* data) const {
        vst1q_u8(data, vector);
    }

    inline NEONVector& operator+=(const NEONVector& rhs) {
        vector = vaddq_u8(vector, rhs.vector);
        return *this;
    }

    inline NEONVector operator==(const NEONVector& rhs) const {
        return NEONVector(vceqq_u8(vector, rhs.vector));
    }

    inline NEONVector operator&(const NEONVector& rhs) const {
        return NEONVector(vandq_u8(vector, rhs.vector));
    }

    inline NEONVector operator|(const NEONVector& rhs) const {
        return NEONVector(vorrq_u8(vector, rhs.vector));
    }

    inline NEONVector operator!() const {
        return NEONVector(vmvnq_u8(vector));
    }

    inline NEONVector andnot(const NEONVector& rhs) const {
        return NEONVector(vbicq_u8(rhs.vector, vector));
    }

    inline uint16_t sum() const {
        #ifdef __aarch64__
            // Don't use `vaddvq_u8` because it can overflow, first add pairwise
            // into 16-bit lanes to avoid overflows
            return vaddvq_u16(vpaddlq_u8(a));
        #else
            uint64x2_t paired = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(a)));
            return vgetq_lane_u64(paired, 0) + vgetq_lane_u64(paired, 1);
        #endif
    }

    inline void clear() {
        vector = vdupq_n_u8(0);
    }
};

namespace statistics {
    void NEONSimilarity::calculateMatrixIdentity() {
        // create a timer for this function
        StartTiming("void NEONSimilarity::calculateMatrixIdentity() ");
        // abort if identity matrix computation was already done
        if (matrixIdentity != nullptr)
            return;
        // Allocate memory for the matrix identity
        matrixIdentity = new float *[alig->originalNumberOfSequences];
        for (int i = 0; i < alig->originalNumberOfSequences; i++) {
            matrixIdentity[i] = new float[alig->originalNumberOfSequences];
        }
        // Run SIMD code with NEON
        simd::calculateMatrixIdentity<NEONVector>(alig, matrixIdentity);
    }

    bool NEONSimilarity::calculateVectors(bool cutByGap) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool NEONSimilarity::calculateVectors(int *gaps) ");

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

        if (!simd::calculateSimilarityVectors<NEONVector>(alig, const_cast<const float**>(matrixIdentity), simMatrix, gaps, MDK)) 
            return false;

        for (int i = 0; i < alig->originalNumberOfSequences; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
        matrixIdentity = nullptr;

        return true;
    }

    void NEONGaps::CalculateVectors() {
        // create a timer for this function
        StartTiming("bool NEONGaps::CalculateVectors(int *gaps) ");
        // calculate gaps in SIMD with NEON
        simd::calculateGapVectors<NEONVector>(alig, gapsInColumn);
        // build histogram and find largest number of gaps
        for (int i = 0; i < alig->originalNumberOfResidues; i++) {
            totalGaps += gapsInColumn[i];
            numColumnsWithGaps[gapsInColumn[i]]++;
            if (gapsInColumn[i] > maxGaps)
                maxGaps = gapsInColumn[i];
        }
    }
}

void NEONCleaner::calculateSeqIdentity() {
  // create a timer for this function
  StartTiming("void NEONCleaner::calculateSeqIdentity(void) ");
  // create identities matrix to store identities scores
  alig->identities = new float*[alig->originalNumberOfSequences];
  for(int i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;
      alig->identities[i] = new float[alig->originalNumberOfSequences];
      alig->identities[i][i] = 0;
  }
  // Run SIMD code with NEON
  simd::calculateSeqIdentity<NEONVector>(alig, alig->identities);
}

bool NEONCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
    // create a timer for this function
    StartTiming("bool NEONCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) ");
    // abort if there is not output vector to write to
    if (spuriousVector == nullptr)
        return false;
    // Run SIMD code with NEON
    return simd::calculateSpuriousVector<NEONVector>(alig, overlap, spuriousVector);
}