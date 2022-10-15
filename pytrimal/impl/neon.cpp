#include <arm_neon.h>
#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

#include "neon.h"
#include "template.h"

class NEONVector {
private:
  uint8x16_t vector;
  inline NEONVector(uint8x16_t vec) : vector(vec) {}

public:
  const static size_t LANES = 16;
  const static size_t SIZE = sizeof(uint8x16_t);

  inline NEONVector() : vector(vdupq_n_u8(0)) {}

  inline static NEONVector duplicate(const uint8_t value) {
    return NEONVector(vdupq_n_u8(value));
  }

  inline static NEONVector load(const uint8_t *data) {
    return NEONVector(vld1q_u8(data));
  }

  inline static NEONVector loadu(const uint8_t *data) {
    return NEONVector(vld1q_u8(data));
  }

  inline void store(uint8_t *data) const { vst1q_u8(data, vector); }

  inline void storeu(uint8_t *data) const { vst1q_u8(data, vector); }

  inline NEONVector &operator+=(const NEONVector &rhs) {
    vector = vaddq_u8(vector, rhs.vector);
    return *this;
  }

  inline NEONVector operator==(const NEONVector &rhs) const {
    return NEONVector(vceqq_u8(vector, rhs.vector));
  }

  inline NEONVector operator&(const NEONVector &rhs) const {
    return NEONVector(vandq_u8(vector, rhs.vector));
  }

  inline NEONVector operator|(const NEONVector &rhs) const {
    return NEONVector(vorrq_u8(vector, rhs.vector));
  }

  inline NEONVector operator!() const { return NEONVector(vmvnq_u8(vector)); }

  inline NEONVector andnot(const NEONVector &rhs) const {
    return NEONVector(vbicq_u8(vector, rhs.vector));
  }

  inline uint16_t sum() const {
#ifdef __aarch64__
    // Don't use `vaddvq_u8` because it can overflow, first add pairwise
    // into 16-bit lanes to avoid overflows
    return vaddvq_u16(vpaddlq_u8(vector));
#else
    uint64x2_t paired = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vector)));
    return vgetq_lane_u64(paired, 0) + vgetq_lane_u64(paired, 1);
#endif
  }

  inline void clear() { vector = vdupq_n_u8(0); }
};

namespace statistics {
void NEONSimilarity::calculateMatrixIdentity() {
  StartTiming("void NEONSimilarity::calculateMatrixIdentity() ");
  simd::calculateMatrixIdentity<NEONVector>(*this);
}

bool NEONSimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool NEONSimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<NEONVector>(*this, cutByGap);
}

void NEONGaps::CalculateVectors() {
  StartTiming("bool NEONGaps::CalculateVectors() ");
  simd::calculateGapVectors<NEONVector>(*this);
}
} // namespace statistics

void NEONCleaner::calculateSeqIdentity() {
  StartTiming("void NEONCleaner::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<NEONVector>(*this);
}

bool NEONCleaner::calculateSpuriousVector(float overlap,
                                          float *spuriousVector) {
  StartTiming("bool NEONCleaner::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<NEONVector>(*this, overlap,
                                                   spuriousVector);
}