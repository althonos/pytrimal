#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Statistics/Similarity.h"
#include "Cleaner.h"

namespace statistics {
    class SSESimilarity: public Similarity {
    private:
        std::vector<char> ascii_vhash;
        std::vector<char> colgap;
        std::string column;
    public:
        SSESimilarity(Alignment* parentAlignment);
        void calculateMatrixIdentity() override;
        bool calculateVectors(bool cutByGap) override;
        bool setSimilarityMatrix(similarityMatrix *sm) override;
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
