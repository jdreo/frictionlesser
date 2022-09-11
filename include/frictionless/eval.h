#pragma once

#include <eo>
#include <frictionless/frictionless.h>
#include <frictionless/transcriptome.h>
#include <frictionless/score.h>

namespace frictionless {

class FullEval : public eoEvalFunc<Signature>
{
    protected:
        FriedmanScore _frs;

    public:

        FullEval(const Transcriptome& ranked, const double alpha=2);

        void operator()(Signature& geneset);
};

} // frictionless
