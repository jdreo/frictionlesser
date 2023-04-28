#include <frictionless/eval.h>

namespace frictionless {

EvalFull::EvalFull(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalFull::operator()(Signature& solution)
{
    CLUTCHLOG(debug, "Full Eval of: " << solution.str());

    ASSERT(_frs.has_transcriptome_cache());
    ASSERT(_frs.has_size_cache());

    _frs.init_signature(solution);

    // // We should have a valid state to swap from.
    // if(not _frs.has_init_signature()) {
    //     CLUTCHLOG(warning, "First initialization of the EvalFull's FriedmanScore");
    //     // _frs.new_signature_size(solution.selected.size());
    //     _frs.init_signature(solution);
    // }

    solution.fitness(Fitness(_frs.score(solution), _frs.swap_cache()));
    // This actually have reset _frs on the given solution.
    ASSERT(not solution.invalid());
}

FriedmanScore& EvalFull::frs()
{
    return _frs;
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
    neighbor.fitness( Fitness(newsol.fitness().score_details(), _full_eval.frs().swap_cache()) );
    ASSERT(not neighbor.invalid());
    // No need to reset _frs, as EvalFull will do it at each call anyway.
}

EvalFull& EvalTest::full_eval()
{
    return _full_eval;
}

#endif


EvalSwap::EvalSwap(FriedmanScore& frs) :
    _frs(frs)
{ }

void EvalSwap::operator()(Signature& solution, moBinaryPartitionSwapNeighbor<Signature> & neighbor)
{
    CLUTCHLOG(xdebug, "Partial eval of: " << solution << " to be moved by: " << neighbor);

    ASSERT(_frs.has_transcriptome_cache());
    ASSERT(_frs.has_size_cache());
    ASSERT(_frs.has_init_signature());

    // // We should have a valid state to swap from.
    // if(not _frs.has_init_signature()) {
    //     CLUTCHLOG(warning, "First initialization of the EvalSwap's FriedmanScore");
    //     // _frs.new_signature_size(solution.selected.size());
    //     _frs.init_signature(solution);
    // }
    ASSERT(not solution.invalid()); // TODO: necessary?
    CLUTCHLOG(xdebug, "Load cache from solution");
    _frs.swap_cache( solution.fitness().cache() );

    CLUTCHLOG(xdebug, "Double check that fitness equals the score...");
    ASSERT(solution.fitness().score_details() == _frs.score(solution));
    CLUTCHLOG(xdebug, "OK");

    // Apply the neighbor move on a temp solution.
    Signature newsol = solution;
    neighbor.move(newsol);
    CLUTCHLOG(xdebug, "Apply to temp neighbor solution: " << newsol);

    // Partial score cache update.
    auto [in,out] = neighbor.get();
    CLUTCHLOG(xdebug, "Prepare swap: " << neighbor);
    _frs.new_swap(in, out);

    // Final score computation and save cache to neighbor.
    Fitness::Type newscore = _frs.score(newsol);
    neighbor.fitness(Fitness(newscore, _frs.swap_cache()));
    CLUTCHLOG(xdebug, "Neighbor solution score: " << neighbor.fitness());
    ASSERT(not neighbor.invalid());

    // Reverse the core state to the one of the origin solution.
    CLUTCHLOG(xdebug, "Reverse swap state: " << neighbor << " to: -" << in << " +" << out);
    // i.e. reload previous cache.
    _frs.swap_cache( solution.fitness().cache() );

    CLUTCHLOG(xdebug, "Previous solution: " << solution << " / score: " << _frs.score(solution));
    CLUTCHLOG(xdebug, "Neighbor solution: " << newsol << " has (neighbor) fitness: " << neighbor.fitness());
    ASSERT(_frs.score(solution) == solution.fitness().score_details());

    ASSERT(not neighbor.invalid());
}

} // frictionless

