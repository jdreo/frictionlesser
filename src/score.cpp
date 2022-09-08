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

    new_transcriptome();

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

void FriedmanScore::clear_cache()
{
    CLUTCHLOG(debug, "Clear cache");
    A     .clear();
    D     .clear();
    B     .clear();
    C     .clear();
    E     .clear();
    F     .clear();
    GG    .clear();
    T_ij  .clear();
    SSR_ij.clear();
    Rc    .clear();

    _has_init = false;
    _cached_signature_size = 0;
    _cached_gene_in = 0;
    _cached_gene_out = 0;
}

void FriedmanScore::new_transcriptome()
{
    CLUTCHLOG(debug, "New transcriptome");
    clear_cache();

    // Precompute constants.
    for(size_t i : _transcriptome.samples() ) { // All samples.
        CLUTCHLOG(xdebug, "    Sample: " << i);
        // Constants depending only on the number of cells
        // in each sample: E, F, GG.
        const double m_i = _transcriptome.cells_nb(i);
        E .push_back( 3 * m_i * std::pow(m_i+1,2) );
        F .push_back( m_i * (m_i+1) );
        GG.push_back( m_i / (m_i - 1) );
        CLUTCHLOG(xdebug, "        E=" << E.back() << ",    F=" << F.back() << ",    GG=" << GG.back());

        // Constants for each cells of each sample: SSR_ij, T_ij.
        // FIXME TODO tests
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : _transcriptome.genes()) { // All genes.

            double sum_ranks_sq = 0; // for SSR_ij.
            // Rank => ties count^3, for T_ij.
            // This should work even if ranks are unordered.
            std::map<double,size_t> ties_cube;

            for(size_t c : _transcriptome.cells(i) ) { // All cells of this sample.
                sum_ranks_sq += std::pow(_transcriptome.rank(c,j),2);
#ifndef NDEBUG
                auto& t = ties_cube[ _transcriptome.rank(c,j) ];
#else
                auto& t = ties_cube.at( _transcriptome.rank(c,j) );
#endif
                t++;

            } // for c in cells
            SSR_j.push_back( sum_ranks_sq );

            // Sum cubic number of ties.
            double sum_ties_cube = 0;
            for([[maybe_unused]] const auto& [r,t] : ties_cube) {
                if(t > 1) { // There were ties for this rank.
                    sum_ties_cube += std::pow(t,3);
            }}
            T_j.push_back( sum_ties_cube );
        } // for j in genes
        SSR_ij.push_back(SSR_j);
        T_ij.push_back(T_j);
        CLUTCHCODE(xdebug,
            std::ostringstream ssr;
            for(double s : SSR_j) {
                ssr << " " << s;
            }
            CLUTCHLOG(xdebug, "        SSR={" << ssr.str() << "}");
            std::ostringstream ts;
            for(double t : T_j) {
                ts << " " << t;
            }
            CLUTCHLOG(xdebug, "        T={" << ts.str() << "}");
        );

    } // for i in samples

}

void FriedmanScore::new_signature_size(const size_t signature_size)
{
    CLUTCHLOG(debug, "New signature size");
    B.clear();
    C.clear();
    ASSERT(signature_size > 0);

    for(size_t i : _transcriptome.samples()) {
        const size_t m_i = _transcriptome.cells_nb(i);
        ASSERT(_transcriptome.cells_nb(i) > 0);
        B.push_back( 3 * std::pow(signature_size,2) * m_i * std::pow(m_i+1, 2) );
        C.push_back(              signature_size    * m_i *         (m_i+1)    );
        CLUTCHLOG(xdebug, "    Sample: " << i);
        CLUTCHLOG(xdebug, "        B=" << B.back() << ",    C=" << C.back());
    } // for i in samples

    _cached_signature_size = signature_size;
    _has_init = false;
}

void FriedmanScore::init_signature(Signature genes)
{
    CLUTCHLOG(debug, "Initialize signature");
    ASSERT(not _has_init);
    ASSERT(A.size() == 0);
    ASSERT(D.size() == 0);
    ASSERT(genes.size() == _transcriptome.genes_nb());

    Rc.clear();
    for(size_t c : _transcriptome.cells()) {
        CLUTCHLOG(xdebug, "    Cell: " << c);
        double sum_r = 0;
        // FIXME use something else than eoBit to improve perf.
        // This should iterate over selected gene indices instead of all genes.
        for(size_t j : _transcriptome.genes()) {
            #ifndef NDEBUG
                const bool selected = genes[j];
            #else
                const bool selected = genes.at(j);
            #endif
            sum_r += selected * _transcriptome.rank(c,j);
            // CLUTCHLOG(xdebug, "         gene=" << j << ", selected=" << selected << ", sum_r=" << sum_r);
        } // for j in genes
        Rc.push_back(sum_r);
        CLUTCHLOG(xdebug, "         Rc=" << Rc.back());
    } // for c in cells

    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOG(xdebug, "    Sample: " << i);
        double sum_Rc2 = 0;
        for(size_t c : _transcriptome.cells(i)) {
            #ifndef NDEBUG
                sum_Rc2 += std::pow(Rc[c], 2);
            #else
                sum_Rc2 += std::pow(Rc.at(c), 2);
            #endif
        } // c in sample cells
        A.push_back(12 * sum_Rc2);

        // HERE
        double sum_t = 0;
        const double m_i = _transcriptome.cells_nb(i);
        for(size_t j : _transcriptome.genes()) {
            #ifndef NDEBUG
                const bool selected = genes[j];
                sum_t += selected * (T_ij[i][j] - m_i);
            #else
                const bool selected = genes.at(j);
                sum_t += selected * (T_ij.at(i).at(j) - m_i);
            #endif
        }
        D.push_back( 1/(m_i-1) * sum_t );

        CLUTCHLOG(xdebug, "        A=" << A.back() << ",    D=" << D.back());

    } // for i in samples

    _has_init = true;
}

void FriedmanScore::new_swap(const size_t gene_in, const size_t gene_out)
{
    CLUTCHLOG(debug, "New swap");
    ASSERT(gene_in != gene_out);

    for(size_t i : _transcriptome.samples()) {
        double sum_Rin = 0;
        double sum_r2 = 0;
        for(size_t c : _transcriptome.cells(i)) {
            // FIXME double check that
            // R_c(G\setminus j_{in}) == R_c(G)-r_{c,j_{in}}
            #ifndef NDEBUG
                sum_Rin +=  _transcriptome.rank(c,gene_out) * (Rc[c] - _transcriptome.rank(c,gene_in))
                          - _transcriptome.rank(c,gene_in ) *  Rc[c];
            #else
                sum_Rin +=  _transcriptome.rank(c,gene_out) * (Rc.at(c) - _transcriptome.rank(c,gene_in))
                          - _transcriptome.rank(c,gene_in ) *  Rc.at(c);
            #endif
            sum_r2 += std::pow(_transcriptome.rank(c,gene_in),2)
                    + std::pow(_transcriptome.rank(c,gene_out),2);
        } // for c in cells
        #ifndef NDEBUG
            const double Ai = A[i];
            A[i] = Ai + 24*sum_Rin + 12*sum_r2;

            const double Di = D[i];
            D[i] = Di + 1/(_transcriptome.cells_nb(i)-1)
                 * ( T_ij[i][gene_out] - T_ij[i][gene_in] );
         #else
            const double Ai = A.at(i);
            A.at(i) = Ai + 24*sum_Rin + 12*sum_r2;

            const double Di = D.at(i);
            D.at(i) = Di + 1/(_transcriptome.cells_nb(i)-1)
                 * ( T_ij.at(i).at(gene_out) - T_ij.at(i).at(gene_in) );
         #endif
    } // for i in samples

    _cached_gene_in = gene_in;
    _cached_gene_out = gene_out;
}


double FriedmanScore::score(Signature genes)
{
    CLUTCHLOG(debug, "Compute score");
    size_t current_signature_size = frictionless::sum(genes);
    ASSERT(current_signature_size > 0);
    if(_cached_signature_size != current_signature_size) {
        new_signature_size(current_signature_size);
    }

    if(not _has_init) {
        init_signature(genes);
    }

    double score = 0;
    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOG(xdebug, "    Sample: " << i);
        double s_hat;
        #ifndef NDEBUG
            double CD = C[i] - D[i];
        #else
            double CD = C.at(i) - D.at(i);
        #endif
        if(CD != 0) {
            #ifndef NDEBUG
                s_hat = ( (A[i]-B[i]) / CD )
            #else
                s_hat = ( (A.at(i)-B.at(i)) / CD )
            #endif
                  / std::pow(_transcriptome.genes_nb(), this->alpha);
        } else {
            s_hat = 0;
        }
        CLUTCHLOG(xdebug, "        s_hat=" << s_hat);
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
