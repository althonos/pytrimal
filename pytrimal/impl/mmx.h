#ifndef _PYTRIMAL_IMPL_MMX
#define _PYTRIMAL_IMPL_MMX

#include "Cleaner.h"
#include "Statistics/Gaps.h"
#include "Statistics/Similarity.h"

namespace statistics {
class MMXSimilarity : public Similarity {
public:
  MMXSimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  void calculateMatrixIdentity() override;
  bool calculateVectors(bool cutByGap) override;
};
class MMXGaps : public Gaps {
public:
  MMXGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  void CalculateVectors() override;
};
} // namespace statistics

class MMXCleaner : public Cleaner {
public:
  MMXCleaner(Alignment *parent) : Cleaner(parent) {}
  void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif
