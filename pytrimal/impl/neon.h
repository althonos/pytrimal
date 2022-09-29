#ifndef _PYTRIMAL_IMPL_NEON
#define _PYTRIMAL_IMPL_NEON

#include "Statistics/Similarity.h"
#include "Statistics/Gaps.h"
#include "Cleaner.h"

namespace statistics {
    class NEONSimilarity: public Similarity {
    public:
        NEONSimilarity(Alignment* parentAlignment): Similarity(parentAlignment) {}
        void calculateMatrixIdentity() override;
        bool calculateVectors(bool cutByGap) override;
    };
    class NEONGaps: public Gaps {
    public:
        NEONGaps(Alignment* parentAlignment): Gaps(parentAlignment) {}
        void CalculateVectors() override;
    };
}

class NEONCleaner: public Cleaner {
public:
    NEONCleaner(Alignment* parent): Cleaner(parent) {}
    void calculateSeqIdentity() override;
    bool calculateSpuriousVector(float overlap, float *spuriousVector);
};

#endif
