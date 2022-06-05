#ifndef _PYTRIMAL_IMPL_SSE
#define _PYTRIMAL_IMPL_SSE

#include "Statistics/Similarity.h"

namespace statistics {

    class SSESimilarity: public Similarity {

    private:
        typedef Similarity super;

    public:
        SSESimilarity(Alignment* parentAlignment);
        void calculateMatrixIdentity() override;
    };

}

#endif
