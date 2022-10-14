#include <climits>
#include <cstdint>
#include <immintrin.h>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

#include "avx.h"
#include "template.h"

class AVXVector {
private:
  __m256i vector;
  inline AVXVector(__m256i vec) : vector(vec) {}

public:
  const static size_t LANES = 32;
  const static size_t SIZE = sizeof(__m256i);

  inline AVXVector() : vector(_mm256_setzero_si256()) {}

  inline static AVXVector duplicate(const uint8_t value) {
    return AVXVector(_mm256_set1_epi8(value));
  }

  inline static AVXVector load(const uint8_t *data) {
    return AVXVector(_mm256_load_si256((const __m256i *)data));
  }

  inline static AVXVector loadu(const uint8_t *data) {
    return AVXVector(_mm256_loadu_si256((const __m256i *)data));
  }

  inline void store(uint8_t *data) const {
    _mm256_store_si256((__m256i *)data, vector);
  }

  inline void storeu(uint8_t *data) const {
    _mm256_storeu_si256((__m256i *)data, vector);
  }

  inline AVXVector &operator+=(const AVXVector &rhs) {
    vector = _mm256_add_epi8(vector, rhs.vector);
    return *this;
  }

  inline AVXVector operator==(const AVXVector &rhs) const {
    return AVXVector(_mm256_cmpeq_epi8(vector, rhs.vector));
  }

  inline AVXVector operator&(const AVXVector &rhs) const {
    return AVXVector(_mm256_and_si256(vector, rhs.vector));
  }

  inline AVXVector operator|(const AVXVector &rhs) const {
    return AVXVector(_mm256_or_si256(vector, rhs.vector));
  }

  inline AVXVector operator!() const {
    return AVXVector(_mm256_andnot_si256(vector, _mm256_set1_epi8(0xFF)));
  }

  inline AVXVector andnot(const AVXVector &rhs) const {
    return AVXVector(_mm256_andnot_si256(rhs.vector, vector));
  }

  inline uint16_t sum() const {
    __m256i vsum = _mm256_sad_epu8(vector, _mm256_setzero_si256());
    return _mm256_extract_epi8(vsum, 0) + _mm256_extract_epi8(vsum, 4) +
           _mm256_extract_epi8(vsum, 8) + _mm256_extract_epi8(vsum, 12);
  }

  inline void clear() { vector = _mm256_setzero_si256(); }
};

namespace statistics {
void AVXSimilarity::calculateMatrixIdentity() {
  StartTiming("void AVXSimilarity::calculateMatrixIdentity() ");
  simd::calculateMatrixIdentity<AVXVector>(*this);
}

bool AVXSimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool AVXSimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<AVXVector>(*this, cutByGap);
}

void AVXGaps::CalculateVectors() {
  StartTiming("bool AVXGaps::CalculateVectors() ");
  simd::calculateGapVectors<AVXVector>(*this);
}
} // namespace statistics

void AVXCleaner::calculateSeqIdentity() {
  StartTiming("void AVXCleaner::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<AVXVector>(*this);
}

bool AVXCleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool AVXCleaner::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<AVXVector>(*this, overlap,
                                                  spuriousVector);
}
