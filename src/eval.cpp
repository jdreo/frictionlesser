#include <frictionless/eval.h>

namespace frictionless {

EvalFull::EvalFull(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalFull::operator()(Signature& geneset)
{
    _frs.new_signature_size(geneset.selected.size());
    _frs.init_signature(geneset);
    geneset.fitness(_frs.score(geneset));
}


EvalSwap::EvalSwap(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalSwap::operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor)
{
    assert(not solution.invalid());

    // Apply the neighbor move on a temp solution.
    Signature newsol = solution;
    neighbor.move(newsol);

    // Partial score cache update.
    auto [in,out] = neighbor.get();
    _frs.new_swap(in, out);

    // Final score computation.
    neighbor.fitness(_frs.score(newsol));
}

} // frictionless

