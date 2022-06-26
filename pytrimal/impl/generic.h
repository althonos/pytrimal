#ifndef _PYTRIMAL_IMPL_GENERIC
#define _PYTRIMAL_IMPL_GENERIC

#include <vector>

#include "Statistics/Similarity.h"

namespace statistics {
    class GenericSimilarity: public Similarity {
    private:
        std::vector<char> column;
        std::vector<char> colgap;
    public:
        GenericSimilarity(Alignment* parentAlignment);
        bool calculateVectors(bool cutByGap) override;
    };
}
#endif
