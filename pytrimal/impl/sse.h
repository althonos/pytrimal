#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Statistics/Similarity.h"
#include "Cleaner.h"

namespace statistics {
    class SSESimilarity: public Similarity {
    public:
        SSESimilarity(Alignment* parentAlignment);
        void calculateMatrixIdentity() override;
    };
}

class SSECleaner: public Cleaner {
private:
    unsigned char* skipResidues;
public:
    SSECleaner(Alignment* parent);
    ~SSECleaner();
    void calculateSeqIdentity() override;
};

#endif
