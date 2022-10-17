#ifndef _PYTRIMAL_IMPL_GENERIC
#define _PYTRIMAL_IMPL_GENERIC

#include <vector>

#include "Cleaner.h"
#include "Statistics/Gaps.h"
#include "Statistics/Similarity.h"

namespace statistics {
class GenericSimilarity : public Similarity {
public:
  GenericSimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  void calculateMatrixIdentity() override;
  bool calculateVectors(bool cutByGap) override;
};
class GenericGaps : public Gaps {
public:
  GenericGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  void CalculateVectors() override;
};
} // namespace statistics

class GenericCleaner : public Cleaner {
public:
  GenericCleaner(Alignment *parent) : Cleaner(parent) {}
  void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif
