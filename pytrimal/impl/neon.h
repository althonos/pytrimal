#ifndef _PYTRIMAL_IMPL_NEON
#define _PYTRIMAL_IMPL_NEON

#include "Cleaner.h"
#include "Statistics/Gaps.h"
#include "Statistics/Similarity.h"

namespace statistics {
class NEONSimilarity : public Similarity {
public:
  NEONSimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  void calculateMatrixIdentity() override;
  bool calculateVectors(bool cutByGap) override;
};
class NEONGaps : public Gaps {
public:
  NEONGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  void CalculateVectors() override;
};
} // namespace statistics

class NEONCleaner : public Cleaner {
public:
  NEONCleaner(Alignment *parent) : Cleaner(parent) {}
  void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif
