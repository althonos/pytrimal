#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Cleaner.h"
#include "Statistics/Gaps.h"
#include "Statistics/Similarity.h"

namespace statistics {
class SSESimilarity : public Similarity {
public:
  SSESimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  void calculateMatrixIdentity() override;
  bool calculateVectors(bool cutByGap) override;
};
class SSEGaps : public Gaps {
public:
  SSEGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  void CalculateVectors() override;
};
} // namespace statistics

class SSECleaner : public Cleaner {
public:
  SSECleaner(Alignment *parent) : Cleaner(parent) {}
  void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif
