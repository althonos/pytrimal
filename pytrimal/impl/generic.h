#ifndef _PYTRIMAL_IMPL_GENERIC
#define _PYTRIMAL_IMPL_GENERIC

#include <vector>

#include "Cleaner.h"
#include "Statistics/Similarity.h"

namespace statistics {
class GenericSimilarity : public Similarity {
public:
  GenericSimilarity(Alignment *parentAlignment);
  bool calculateVectors(bool cutByGap) override;
};
} // namespace statistics

class GenericCleaner : public Cleaner {
private:
  // temporary counter for `calculateSpuriousVector`
  uint32_t *hits;

public:
  GenericCleaner(Alignment *parent);
  ~GenericCleaner();
  // void calculateSeqIdentity() override;
  bool calculateSpuriousVector(float overlap, float *spuriousVector);
};

#endif
