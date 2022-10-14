#include <climits>
#include <cstdint>
#include <mmintrin.h>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

#include "mmx.h"
#include "template.h"

class MMXVector {
private:
  __m64 vector;
  inline MMXVector(__m64 vec) : vector(vec) {}

public:
  const static size_t LANES = 8;
  const static size_t SIZE = sizeof(__m64);

  inline MMXVector() : vector(_mm_setzero_si64()) {}

  inline static MMXVector duplicate(const uint8_t value) {
    return MMXVector(_mm_set1_pi8(value));
  }

  inline static MMXVector load(const uint8_t* data) {
    return MMXVector(*((const __m64*) data));
  }

  inline static MMXVector loadu(const uint8_t* data) {
    return MMXVector(*((const __m64*) data));
  }

  inline void store(uint8_t *data) const {
    *((__m64*) data) = vector;
  }

  inline void storeu(uint8_t *data) const {
    *((__m64*) data) = vector;
  }

  inline MMXVector &operator+=(const MMXVector &rhs) {
    vector = _mm_add_pi8(vector, rhs.vector);
    return *this;
  }

  inline MMXVector operator==(const MMXVector &rhs) const {
    return MMXVector(_mm_cmpeq_pi8(vector, rhs.vector));
  }

  inline MMXVector operator&(const MMXVector &rhs) const {
    return MMXVector(_mm_and_si64(vector, rhs.vector));
  }

  inline MMXVector operator|(const MMXVector &rhs) const {
    return MMXVector(_mm_or_si64(vector, rhs.vector));
  }

  inline MMXVector operator!() const {
    return MMXVector(_mm_andnot_si64(vector, _mm_set1_pi8(0xFF)));
  }

  inline MMXVector andnot(const MMXVector &rhs) const {
    return MMXVector(_mm_andnot_si64(rhs.vector, vector));
  }

  inline uint16_t sum() const {
    __m64    hi      = _mm_unpackhi_pi16(_mm_setzero_si64(), vector);
    __m64    lo      = _mm_unpacklo_pi16(_mm_setzero_si64(), vector);
    __m64    partial = _mm_add_pi16(hi, lo);
    uint64_t value   = _mm_cvtm64_si64(partial);
    return (value & 0xFFFF) + (value >> 16);
  }

  inline void clear() { vector = _mm_setzero_si64(); }
};

namespace statistics {
void MMXSimilarity::calculateMatrixIdentity() {
  StartTiming("void MMXSimilarity::calculateMatrixIdentity() ");
  simd::calculateMatrixIdentity<MMXVector>(*this);
}

bool MMXSimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool MMXSimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<MMXVector>(*this, cutByGap);
}

void MMXGaps::CalculateVectors() {
  StartTiming("bool MMXGaps::CalculateVectors() ");
  simd::calculateGapVectors<MMXVector>(*this);
}
} // namespace statistics

void MMXCleaner::calculateSeqIdentity() {
  StartTiming("void MMXCleaner::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<MMXVector>(*this);
}

bool MMXCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool MMXCleaner::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<MMXVector>(*this, overlap,
                                                  spuriousVector);
}
