#include <frictionless/eval.h>

namespace frictionless {

EvalFull::EvalFull(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalFull::operator()(Signature& solution)
{
    CLUTCHLOG(debug, "Full Eval of: " << solution);
    _frs.new_signature_size(solution.selected.size());
    _frs.init_signature(solution);
    solution.fitness(_frs.score(solution));
    // This actually have reset _frs on the given solution.
    ASSERT(not solution.invalid());
}


#ifndef NDEBUG
EvalTest::EvalTest(EvalFull& full_eval) :
    _full_eval(full_eval)
{ }

void EvalTest::operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor)
{
    CLUTCHLOG(debug, "Full Eval of: " << solution << " moved as " << neighbor);
    // Apply the neighbor move on a temp solution.
    Signature newsol = solution;
    neighbor.move(newsol);
    // Full evaluation.
    _full_eval(newsol);
    neighbor.fitness(newsol.fitness());
    ASSERT(not neighbor.invalid());
    // No need to reset _frs, as EvalFull will do it at each call anyway.
}
#endif


EvalSwap::EvalSwap(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalSwap::operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor)
{
    CLUTCHLOG(debug, "Partial eval of: " << solution << " moved as " << neighbor);

    // We should have a valid state to swap from.
    if(not _frs.has_init_signature()) {
        CLUTCHLOG(xdebug, "First initialization of the EvalSwap's FriedmanScore");
        _frs.new_signature_size(solution.selected.size());
        _frs.init_signature(solution);
    }
    // ASSERT(not solution.invalid());
    // ASSERT(solution.fitness() == _frs.score(solution));

    // Save the partial state of the given solution.
    std::vector<double> A = _frs.A;
    std::vector<double> D = _frs.D;
    std::vector<double> R = _frs.R;

    // Apply the neighbor move on a temp solution.
    Signature newsol = solution;
    neighbor.move(newsol);
    CLUTCHLOG(debug, "    That is solution: " << newsol);

    // Partial score cache update.
    auto [in,out] = neighbor.get();
    CLUTCHLOG(debug, "Prepare swap: " << neighbor);
    _frs.new_swap(in, out);

    // Final score computation.
    neighbor.fitness(_frs.score(newsol));
    CLUTCHLOG(debug, "Score: " << neighbor.fitness());
    ASSERT(not neighbor.invalid());

    // Reverse the core state to the one of the origin solution.
    CLUTCHLOG(debug, "Reverse swap");
    _frs.A = std::move(A);
    _frs.D = std::move(D);
    _frs.R = std::move(R);
    _frs._cached_gene_in  = std::move(out);
    _frs._cached_gene_out = std::move(in);
    CLUTCHLOG(debug, "Solution: " << solution << " has fitness: " << solution.fitness() << " / score: " << _frs.score(solution));
    ASSERT(_frs.score(solution) == solution.fitness());
}

} // frictionless

