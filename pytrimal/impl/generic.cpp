#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

#include "generic.h"
#include "template.h"

class GenericVector {
private:
  uint8_t vector;
  inline GenericVector(const uint8_t value) : vector(value) {}

public:
  const static size_t LANES = 1;
  const static size_t SIZE = sizeof(uint8_t);

  inline GenericVector() : vector(0) {}

  inline static GenericVector duplicate(const uint8_t value) {
    return GenericVector(value);
  }

  inline static GenericVector load(const uint8_t *data) {
    return GenericVector(*data);
  }

  inline static GenericVector loadu(const uint8_t *data) {
    return GenericVector(*data);
  }

  inline void store(uint8_t *data) const {
    *data = vector;
  }

  inline void storeu(uint8_t *data) const {
    *data = vector;
  }

  inline GenericVector &operator+=(const GenericVector &rhs) {
    vector += rhs.vector;
    return *this;
  }

  inline GenericVector operator==(const GenericVector &rhs) const {
    return GenericVector((vector == rhs.vector) ? 0xFF : 0);
  }

  inline GenericVector operator&(const GenericVector &rhs) const {
    return GenericVector(vector & rhs.vector);
  }

  inline GenericVector operator|(const GenericVector &rhs) const {
    return GenericVector(vector | rhs.vector);
  }

  inline GenericVector operator!() const {
    return GenericVector(~vector);
  }

  inline GenericVector andnot(const GenericVector &rhs) const {
    return GenericVector(vector & (~rhs.vector));
  }

  inline uint16_t sum() const {
    return vector;
  }

  inline void clear() {
    vector = 0;
  }
};

namespace statistics {
void GenericSimilarity::calculateMatrixIdentity() {
  StartTiming("void GenericSimilarity::calculateMatrixIdentity() ");
  simd::calculateMatrixIdentity<GenericVector>(*this);
}

bool GenericSimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool GenericSimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<GenericVector>(*this, cutByGap);
}

void GenericGaps::CalculateVectors() {
  StartTiming("bool GenericGaps::CalculateVectors() ");
  simd::calculateGapVectors<GenericVector>(*this);
}
} // namespace statistics

void GenericCleaner::calculateSeqIdentity() {
  StartTiming("void GenericCleaner::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<GenericVector>(*this);
}

bool GenericCleaner::calculateSpuriousVector(float overlap,
                                          float *spuriousVector) {
  StartTiming("bool GenericCleaner::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<GenericVector>(*this, overlap,
                                                   spuriousVector);
}