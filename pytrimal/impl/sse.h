#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Statistics/Similarity.h"
#include "Cleaner.h"

namespace statistics {
    class SSESimilarity: public Similarity {
    public:
        SSESimilarity(Alignment* parentAlignment);
        void calculateMatrixIdentity() override;
        bool calculateVectors(bool cutByGap) override;
    };
    class SSEGaps: public Gaps {
    public:
        SSEGaps(Alignment* parentAlignment);
        void CalculateVectors() override;
    };
}

class SSECleaner: public Cleaner {
private:
    // SIMD index for which residues can be skipped in `skipResidues`
    unsigned char* skipResidues;
    unsigned char* skipResidues_unaligned;
    // temporary counters for `calculateSpuriousVector`
    uint8_t* hits_u8;
    uint8_t* hits_u8_unaligned;
    uint32_t* hits;
    uint32_t* hits_unaligned;
public:
    SSECleaner(Alignment* parent);
    ~SSECleaner();
    void calculateSeqIdentity() override;
    bool calculateSpuriousVector(float overlap, float *spuriousVector);
};

#endif
