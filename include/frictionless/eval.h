#pragma once

#include <eo>
#include <frictionless/frictionless.h>
#include <frictionless/transcriptome.h>
#include <frictionless/score.h>

namespace frictionless {

/** Full evaluator of a @ref Signature.
 *
 * Do not perform any partial update.
 */
class FullEval : public eoEvalFunc<Signature>
{
    protected:
        /** Score function and its cache. */
        FriedmanScore _frs;

    public:
        /** Constructor.
         *
         * @param ranked A ranked expression matrix.
         * @param alpha The alpha parameter of the score function.
         */
        FullEval(const Transcriptome& ranked, const double alpha=2);

        /** Compute the score function and set the result as a fitness of the given signature. */
        void operator()(Signature& geneset);
};

} // frictionless
