#ifndef _PYTRIMAL_IMPL_AVX
#define _PYTRIMAL_IMPL_AVX

#include "Cleaner.h"
#include "Statistics/Gaps.h"
#include "Statistics/Similarity.h"

namespace statistics {
class AVXSimilarity : public Similarity {
public:
  AVXSimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  void calculateMatrixIdentity() override;
  bool calculateVectors(bool cutByGap) override;
};
class AVXGaps : public Gaps {
public:
  AVXGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  void CalculateVectors() override;
};
} // namespace statistics

class AVXCleaner : public Cleaner {
public:
  AVXCleaner(Alignment *parent) : Cleaner(parent) {}
  void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif