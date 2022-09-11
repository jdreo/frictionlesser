#include <frictionless/eval.h>

namespace frictionless {

FullEval::FullEval(const Transcriptome& ranked, const double alpha) :
    _frs(ranked, alpha)
{ }

void FullEval::operator()(Signature& geneset)
{
    _frs.new_signature_size(sum(geneset));
    _frs.init_signature(geneset);
    geneset.fitness(_frs.score(geneset));
}

} // frictionless

