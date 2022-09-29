#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Statistics/Similarity.h"
#include "Statistics/Gaps.h"
#include "Cleaner.h"

namespace statistics {
    class SSESimilarity: public Similarity {
    public:
        SSESimilarity(Alignment* parentAlignment): Similarity(parentAlignment) {}
        void calculateMatrixIdentity() override;
        bool calculateVectors(bool cutByGap) override;
    };
    class SSEGaps: public Gaps {
    public:
        SSEGaps(Alignment* parentAlignment): Gaps(parentAlignment) {}
        void CalculateVectors() override;
    };
}

class SSECleaner: public Cleaner {
public:
    SSECleaner(Alignment* parent): Cleaner(parent) {}
    void calculateSeqIdentity() override;
    bool calculateSpuriousVector(float overlap, float *spuriousVector);
};

#endif
