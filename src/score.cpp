#include <cmath>
#include <vector>

#include "frictionless/frictionless.h"
#include "frictionless/score.h"
#include "R/pgamma.h"

namespace frictionless {

FriedmanScore::FriedmanScore( const Transcriptome& rt, const double a) :
    _transcriptome(rt),
    alpha(a),
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
    T  .reserve(samples_nb);
    SSR.reserve(samples_nb);

    Rc.reserve(_transcriptome.cells_nb());

    new_transcriptome();

#ifndef NDEBUG
    // Basic size checks.
    ASSERT(E .size() == samples_nb);
    ASSERT(F .size() == samples_nb);
    ASSERT(GG.size() == samples_nb);
    ASSERT(SSR.size() == samples_nb);
    ASSERT(  T.size() == samples_nb);
    for(auto row : SSR) {
        ASSERT(row.size() == _transcriptome.genes_nb()); }
    for(auto row : T) {
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
    T  .clear();
    SSR.clear();
    Rc    .clear();

    _cached_signature_size = 0;
    _cached_gene_in = 0;
    _cached_gene_out = 0;
}

/******************************************
 * E, F, GG, SSR, T
 ******************************************/
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

        // Constants for each cells of each sample: SSR, T.
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : _transcriptome.genes()) { // All genes.

            double sum_ranks_sq = 0; // for SSR.
            // Rank => ties count^3, for T.
            // This should work even if ranks are unordered.
            std::map<double,size_t> ties; // For T.
            for(size_t c : _transcriptome.cells(i) ) { // All cells of this sample.
                const double rk = _transcriptome.rank(c,j);
                sum_ranks_sq += std::pow(rk,2); // For SSR.
                if(ties.contains(rk)) { // For T.
                    ties[rk]++;
                } else {
                    ties[rk] = 1;
                }
            } // for c in cells
            SSR_j.push_back( sum_ranks_sq );

            // Sum cubic number of ties.
            double sum_ties_cube = 0;
            for([[maybe_unused]] const auto& [r,t] : ties) {
                sum_ties_cube += std::pow(t,3);
            }
            T_j.push_back( sum_ties_cube );

        } // for j in genes
        SSR.push_back(SSR_j);
        T.push_back(T_j);
        CLUTCHCODE(xdebug,
            std::ostringstream ssr;
            for(double s : SSR_j) { ssr << " " << s; }
            CLUTCHLOG(xdebug, "        SSR={" << ssr.str() << "}");
            std::ostringstream ts;
            for(double t : T_j  ) {  ts << " " << t; }
            CLUTCHLOG(xdebug, "        T={" << ts.str() << "}");
        );

    } // for i in samples

}

/******************************************
 * B, C
 ******************************************/
void FriedmanScore::new_signature_size(const size_t signature_size)
{
    CLUTCHLOG(debug, "New signature size: " << signature_size);
    B.clear();
    C.clear();
    ASSERT(signature_size > 0);

    for(size_t i : _transcriptome.samples()) {
        ASSERT(_transcriptome.cells_nb(i) > 0);
        #ifndef NDEBUG
            B.push_back( std::pow(signature_size,2) * E[i] );
            C.push_back(          signature_size    * F[i]);
        #else
            B.push_back( std::pow(signature_size,2) * E.at(i) );
            C.push_back(          signature_size    * F.at(i));
        #endif
        CLUTCHLOG(xdebug, "    Sample: " << i);
        CLUTCHLOG(xdebug, "        B=" << B.back() << ",    C=" << C.back());
    } // for i in samples

    _cached_signature_size = signature_size;
}

/******************************************
 * A, D
 ******************************************/
void FriedmanScore::new_swap(const size_t gene_in, const size_t gene_out)
{
    CLUTCHLOG(debug, "New swap: in=" << gene_in << ", out=" << gene_out);
    ASSERT(gene_in != gene_out);
    ASSERT(not (gene_in == _cached_gene_in and gene_out == _cached_gene_out));

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
                 * ( T[i][gene_out] - T[i][gene_in] );
         #else
            const double Ai = A.at(i);
            A.at(i) = Ai + 24*sum_Rin + 12*sum_r2;

            const double Di = D.at(i);
            D.at(i) = Di + 1/(_transcriptome.cells_nb(i)-1)
                 * ( T.at(i).at(gene_out) - T.at(i).at(gene_in) );
         #endif
    } // for i in samples

    _cached_gene_in = gene_in;
    _cached_gene_out = gene_out;
}

/******************************************
 * Rc, A, D
 ******************************************/
void FriedmanScore::init_signature(Signature genes)
{
    CLUTCHLOG(debug, "Initialize signature: " << genes);
    ASSERT(A.size() == 0);
    ASSERT(D.size() == 0);
    ASSERT(genes.size() == _transcriptome.genes_nb());
    ASSERT(_cached_signature_size == sum(genes));

    Rc.clear();
    for(size_t c : _transcriptome.cells()) {
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
        CLUTCHLOG(xdebug, "         R_" << c << "=" << Rc.back());
    } // for c in cells

    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOG(xdebug, "    Sample: " << i);
        double sum_Rc2 = 0;
        for(size_t c : _transcriptome.cells(i)) {
            #ifndef NDEBUG
                sum_Rc2 += 12*std::pow(Rc[c], 2);
            #else
                sum_Rc2 += 12*std::pow(Rc.at(c), 2);
            #endif
        } // c in sample cells
        // A.push_back(12 * sum_Rc2);
        A.push_back(sum_Rc2);

        // Paper version:
        // double sum_t = 0;
        // const double m_i = _transcriptome.cells_nb(i);
        // for(size_t j : _transcriptome.genes()) {
        //     #ifndef NDEBUG
        //         const bool selected = genes[j];
        //         sum_t += selected * (T[i][j] - m_i);
        //     #else
        //         const bool selected = genes.at(j);
        //         sum_t += selected * (T.at(i).at(j) - m_i);
        //     #endif
        // }
        // D.push_back( 1/(m_i-1) * sum_t );

        // Zakiev code:
        double sum_t = 0;
        for(size_t j : _transcriptome.genes()) {
            const bool selected = genes[j];
            sum_t += selected * T[i][j];
        }
        const double m_i = _transcriptome.cells_nb(i);
        D.push_back( (1/(m_i-1)) * sum_t - GG[i] * _cached_signature_size );
        // HERE

        CLUTCHLOG(xdebug, "        A=" << A.back() << ",    D=" << D.back());

    } // for i in samples
}


double FriedmanScore::score(Signature genes)
{
    CLUTCHLOG(debug, "Compute score of: " << genes);
    size_t current_signature_size = frictionless::sum(genes);
    ASSERT(current_signature_size > 0);
    if(_cached_signature_size != current_signature_size) {
        new_signature_size(current_signature_size);
    }

    double score = 0;
    for(size_t i : _transcriptome.samples()) {
        CLUTCHLOG(xdebug, "    Sample: " << i);
        double s_hat;
        #ifndef NDEBUG
            const double CD = C[i] - D[i];
        #else
            const double CD = C.at(i) - D.at(i);
        #endif
        if(CD != 0) {
            #ifndef NDEBUG
                s_hat = ( (A[i]-B[i]) / CD )
            #else
                s_hat = ( (A.at(i)-B.at(i)) / CD )
            #endif
                  / std::pow(current_signature_size, this->alpha);
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

    ASSERT(not std::isnan(score));
    ASSERT(score >= 0);
    ASSERT(not std::isinf(score));
    return score;
}

double FriedmanScore::sqrt_logchisq(const double s, const double m) const
{
    return std::sqrt(-1 * R::pgamma(/*x*/s, /*alph*/m/2, /*scale*/2, /*lower_tail*/0, /*log_p*/1));
}

} // frictionless
