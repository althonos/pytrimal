#include <climits>
#include <cstdint>
#include <emmintrin.h>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

#include "sse.h"
#include "template.h"

class SSEVector {
private:
  __m128i vector;
  inline SSEVector(__m128i vec) : vector(vec) {}

public:
  const static size_t LANES = 16;
  const static size_t SIZE = sizeof(__m128i);

  inline SSEVector() : vector(_mm_setzero_si128()) {}

  inline static SSEVector duplicate(const uint8_t value) {
    return SSEVector(_mm_set1_epi8(value));
  }

  inline static SSEVector load(const uint8_t *data) {
    return SSEVector(_mm_load_si128((const __m128i *)data));
  }

  inline static SSEVector loadu(const uint8_t *data) {
    return SSEVector(_mm_loadu_si128((const __m128i *)data));
  }

  inline void store(uint8_t *data) const {
    _mm_store_si128((__m128i *)data, vector);
  }

  inline void storeu(uint8_t *data) const {
    _mm_storeu_si128((__m128i *)data, vector);
  }

  inline SSEVector &operator+=(const SSEVector &rhs) {
    vector = _mm_add_epi8(vector, rhs.vector);
    return *this;
  }

  inline SSEVector operator==(const SSEVector &rhs) const {
    return SSEVector(_mm_cmpeq_epi8(vector, rhs.vector));
  }

  inline SSEVector operator&(const SSEVector &rhs) const {
    return SSEVector(_mm_and_si128(vector, rhs.vector));
  }

  inline SSEVector operator|(const SSEVector &rhs) const {
    return SSEVector(_mm_or_si128(vector, rhs.vector));
  }

  inline SSEVector operator!() const {
    return SSEVector(_mm_andnot_si128(vector, _mm_set1_epi8(0xFF)));
  }

  inline SSEVector andnot(const SSEVector &rhs) const {
    return SSEVector(_mm_andnot_si128(rhs.vector, vector));
  }

  inline uint16_t sum() const {
    __m128i vsum = _mm_sad_epu8(vector, _mm_setzero_si128());
    return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
  }

  inline void clear() { vector = _mm_setzero_si128(); }
};

namespace statistics {
void SSESimilarity::calculateMatrixIdentity() {
  StartTiming("void SSESimilarity::calculateMatrixIdentity() ");
  simd::calculateMatrixIdentity<SSEVector>(*this);
}

bool SSESimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool SSESimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<SSEVector>(*this, cutByGap);
}

void SSEGaps::CalculateVectors() {
  StartTiming("bool SSEGaps::CalculateVectors() ");
  simd::calculateGapVectors<SSEVector>(*this);
}
} // namespace statistics

void SSECleaner::calculateSeqIdentity() {
  StartTiming("void SSECleaner::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<SSEVector>(*this);
}

bool SSECleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool SSECleaner::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<SSEVector>(*this, overlap,
                                                  spuriousVector);
}
