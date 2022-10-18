#pragma once

#include <eo>
#include <frictionless/frictionless.h>
#include <frictionless/transcriptome.h>
#include <frictionless/score.h>
#include <frictionless/moBinaryPartitionSwapNeighbor.h>

namespace frictionless {

/** Full evaluator of a @ref Signature.
 *
 * Do not perform any partial update.
 */
class EvalFull : public eoEvalFunc<Signature>
{
    protected:
        /** Score function and its cache. */
        FriedmanScore& _frs;

    public:
        /** Constructor.
         *
         * @param ranked A ranked expression matrix.
         * @param alpha The alpha parameter of the score function.
         */
        EvalFull(FriedmanScore& frs);

        /** Compute the score function and set the result as a fitness of the given signature. */
        void operator()(Signature& geneset);
};

#ifndef NDEBUG
/** Partial evaluator that actually perform a full evaluation.
 *
 * This is to test that a partial evaluation ends on the same value than a full evaluation.
 */
class EvalTest : public moEval<moBinaryPartitionSwapNeighbor<Signature>>
{
    protected:
        EvalFull& _full_eval;

    public:
        /** Constructor.
         *
         * @param full_eval The full evaluator.
         */
        EvalTest(EvalFull& full_eval);

        /** Incremental evaluation.
         *
         * @param solution the solution on which to apply the move.
         * @param neighbor the move to consider.
         */
        void operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor); 
};
#endif

/** Partial evaluator of a neighbor that is one swap away.
 */
class EvalSwap : public moEval<moBinaryPartitionSwapNeighbor<Signature>>
{
    protected:
        /** Score function and its cache. */
        FriedmanScore& _frs;

    public:
        /** Constructor.
         *
         * @param ranked A ranked expression matrix.
         * @param alpha The alpha parameter of the score function.
         */
        EvalSwap(FriedmanScore& frs);

        /** Incremental evaluation.
         *
         * @param solution the solution on which to apply the move.
         * @param neighbor the move to consider.
         */
        void operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor); 
};

} // frictionless
