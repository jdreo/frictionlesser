#include <cmath>
#include <vector>

#include "frictionless/frictionless.h"
#include "frictionless/score.h"
#include "R/pgamma.h"

namespace frictionless {

FriedmanScore::FriedmanScore( const Transcriptome& rt, const double a) :
    _transcriptome(rt),
    alpha(a),
    _has_init(false),
    _cached_signature_size(0)
{
    ASSERT(_transcriptome.ranks( ).size() > 0);
    ASSERT(_transcriptome.ranks(0).size() > 0);
    ASSERT(_transcriptome.genes_nb() > 0);
    ASSERT(_transcriptome.samples_nb() > 0);
    ASSERT(_transcriptome.cells_nb(0) > 0);
    ASSERT(alpha > 0);

    const size_t samples_nb = _transcriptome.samples_nb();
    A     .reserve(samples_nb);
    D     .reserve(samples_nb);
    B     .reserve(samples_nb);
    C     .reserve(samples_nb);
    E     .reserve(samples_nb);
    F     .reserve(samples_nb);
    GG    .reserve(samples_nb);
    T_ij  .reserve(samples_nb);
    SSR_ij.reserve(samples_nb);

    Rc.reserve(_transcriptome.cells_nb());

    cache_transcriptome();

#ifndef NDEBUG
    // Basic size checks.
    ASSERT(E .size() == samples_nb);
    ASSERT(F .size() == samples_nb);
    ASSERT(GG.size() == samples_nb);
    ASSERT(SSR_ij.size() == samples_nb);
    ASSERT(  T_ij.size() == samples_nb);
    for(auto row : SSR_ij) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
    for(auto row : T_ij) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }

#endif
}

void FriedmanScore::cache_transcriptome()
{
    E.clear();
    F.clear();
    GG.clear();
    SSR_ij.clear();
    T_ij.clear();

    // Precompute constants.
    for(size_t i : _transcriptome.samples() ) { // All samples.
        // Constants depending only on the number of cells
        // in each sample: E, F, GG.
        const double m_i = _transcriptome.cells_nb(i);
        E .push_back( 3 * m_i * std::pow(m_i+1,2) );
        F .push_back( m_i * (m_i+1) );
        GG.push_back( m_i / (m_i - 1) );

        // Constants for each cells of each sample: SSR_ij, T_ij.
        // FIXME TODO tests
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : _transcriptome.genes()) { // All genes.

            double sum_ranks_sq = 0; // for SSR_ij.
            // Rank => ties count^3, for T_ij.
            // This should work even if ranks are unordered.
            std::map<double,std::vector<double>> ties_cube;

            for(size_t c : _transcriptome.cells(i) ) { // All cells of this sample.
                sum_ranks_sq += std::pow(_transcriptome.rank(c,j),2);
#ifndef NDEBUG
                auto& t = ties_cube[ _transcriptome.rank(c,j) ];
#else
                auto& t = ties_cube.at( _transcriptome.rank(c,j) );
#endif
                t.push_back( std::pow(
                    1 + t.size(), 3
                ));

            } // for c in cells
            SSR_j.push_back( sum_ranks_sq );

            // Sum cubic number of ties.
            double sum_ties_cube = 0;
            for([[maybe_unused]] const auto& [r,t] : ties_cube) {
                if(t.size() > 1) { // There were ties for this rank.
                    sum_ties_cube += sum(t);
            }}
            T_j.push_back( sum_ties_cube );
        } // for j in genes
        SSR_ij.push_back(SSR_j);
        T_ij.push_back(T_j);

    } // for i in samples
}

void FriedmanScore::cache_signature_size(const size_t signature_size)
{
    B.clear();
    C.clear();

    for(size_t i : _transcriptome.samples()) {
        B.push_back( 3 * std::pow(signature_size,2) * _transcriptome.cells_nb(i) * std::pow(_transcriptome.cells_nb(i)+1, 2) );
        C.push_back(              signature_size    * _transcriptome.cells_nb(i) *          _transcriptome.cells_nb(i)+1     );
    } // for i in samples

    _cached_signature_size = signature_size;
}

void FriedmanScore::init_signature(Signature genes)
{
    ASSERT(_has_init = false);
    ASSERT(A.size() == 0);
    ASSERT(D.size() == 0);

    Rc.clear();
    for(size_t c : _transcriptome.cells()) {
        double sum_r = 0;
        for(size_t j : _transcriptome.genes()) {
            for(bool selected : genes) {
                sum_r += selected * _transcriptome.rank(c,j);
            }
        } // for j in genes
        Rc.push_back(sum_r);
    } // for c in cells

    for(size_t i : _transcriptome.samples()) {
        double sum_Rc2 = 0;
        for(size_t c : _transcriptome.cells(i)) {
            // FIXME use something else than eoBit to improve perf.
            // This should iterate over selected gene indices instead of all genes.
            for(size_t j=0; j < genes.size(); ++j) {
                const bool selected = genes[j];
                sum_Rc2 += selected * _transcriptome.rank(c,j);
            }
        } // c in sample cells
        A.push_back(12 * sum_Rc2);

        double sum_t = 0;
        for(size_t j=0; j < genes.size(); ++j) {
            const bool selected = genes[j];
            sum_t += selected * (T_ij[i][j] - _transcriptome.cells_nb(i));
        }
        D.push_back( 1/(_transcriptome.cells_nb(i)-1) * sum_t );

    } // for i in samples

    _has_init = true;
}

void FriedmanScore::cache_swap(const size_t gene_in, const size_t gene_out)
{
    double sum_Rin = 0;
    double sum_r2 = 0;
    for(size_t i : _transcriptome.samples()) {
        for(size_t c : _transcriptome.cells(i)) {
            // FIXME double check that
            // R_c(G\setminus j_{in}) == R_c(G)-r_{c,j_{in}}
            sum_Rin +=  _transcriptome.rank(c,gene_out) * (Rc[c] - _transcriptome.rank(c,gene_in))
                      - _transcriptome.rank(c,gene_in ) *  Rc[c];
            sum_r2 += std::pow(_transcriptome.rank(c,gene_in),2)
                    + std::pow(_transcriptome.rank(c,gene_out),2);
        }
        const double Ai = A[i];
        A[i] = Ai + 24*sum_Rin + 12*sum_r2;

        const double Di = D[i];
        D[i] = Di + 1/(_transcriptome.cells_nb(i)-1)
             * ( T_ij[i][gene_out] - T_ij[i][gene_in] );
    }

    _cached_gene_in = gene_in;
    _cached_gene_out = gene_out;
}


double FriedmanScore::score(Signature genes)
{
    size_t current_signature_size = frictionless::sum(genes);
    if(_cached_signature_size != current_signature_size) {
        cache_signature_size(current_signature_size);
    }

    if(not _has_init) {
        init_signature(genes);
    }

    double score = 0;
    for(size_t i : _transcriptome.samples()) {
        double s_hat;
        double CD = C[i] - D[i];
        if(CD != 0) {
            s_hat = ( (A[i]-B[i]) / CD )
                  / std::pow(_transcriptome.genes_nb(), this->alpha);
        } else {
            s_hat = 0;
        }
        ASSERT(not std::isnan(s_hat));
        ASSERT(s_hat >= 0);
        ASSERT(not std::isinf(s_hat));

        double logpval = sqrt_logchisq(s_hat, _transcriptome.cells_nb(i)-1);

        score += logpval;
    }

    return score;
}

double FriedmanScore::sqrt_logchisq(const double s, const double m) const
{
    return std::sqrt(-1 * R::pgamma(/*x*/s, /*alph*/m/2, /*scale*/2, /*lower_tail*/0, /*log_p*/1));
}

} // frictionless
