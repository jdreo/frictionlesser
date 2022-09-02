#include <cmath>
#include <vector>

#include "frictionless/frictionless.h"
#include "frictionless/score.h"

namespace frictionless {

FriedmanScore::FriedmanScore( const Transcriptome& rt, const double a)
    : _transcriptome(rt), alpha(a)
{
    ASSERT(rt.ranks( ).size() > 0);
    ASSERT(rt.ranks(0).size() > 0);
    ASSERT(rt.genes_nb() > 0);
    ASSERT(rt.samples_nb() > 0);
    ASSERT(rt.cells_nb(0) > 0);
    ASSERT(alpha > 0);

    const size_t samples_nb = rt.samples_nb();
    A     .reserve(samples_nb);
    D     .reserve(samples_nb);
    B     .reserve(samples_nb);
    C     .reserve(samples_nb);
    E     .reserve(samples_nb);
    F     .reserve(samples_nb);
    GG    .reserve(samples_nb);
    T_ij  .reserve(samples_nb);
    SSR_ij.reserve(samples_nb);

    // Precompute constants.
    for(size_t i : rt.samples() ) { // All samples.
        // Constants depending only on the number of cells
        // in each sample.
        const double m_i = rt.cells_nb(i);
        E .push_back( 3 * m_i * std::pow(m_i+1,2) );
        F .push_back( m_i * (m_i+1) );
        GG.push_back( m_i / (m_i - 1) );

        // Constants for each cells of each sample.
        // FIXME TODO tests
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : rt.genes()) { // All genes.

            double sum_ranks_sq = 0; // for SSR_ij.
            // Rank => ties count^3, for T_ij.
            // This should work even if ranks are unordered.
            std::map<double,std::vector<double>> ties_cube;

            for(size_t c : rt.cells(i) ) { // All cells of this sample.
                sum_ranks_sq += std::pow(rt.rank(c,j),2);
#ifndef NDEBUG
                auto& t = ties_cube[ rt.rank(c,j) ];
#else
                auto& t = ties_cube.at( rt.rank(c,j) );
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
#ifndef NDEBUG
            T_j.push_back( sum_ties_cube );
#else
            T_j.push_back( sum_ties_cube );
#endif
        } // for j in genes
        SSR_ij.push_back(SSR_j);
        T_ij.push_back(T_j);
    } // for i in samples

#ifndef NDEBUG
    // Basic size checks.
    ASSERT(E .size() == samples_nb);
    ASSERT(F .size() == samples_nb);
    ASSERT(GG.size() == samples_nb);
    ASSERT(SSR_ij.size() == samples_nb);
    ASSERT(  T_ij.size() == samples_nb);
    for(auto row : SSR_ij) {
        ASSERT(row.size() == rt.genes_nb()); }
    for(auto row : T_ij) {
        ASSERT(row.size() == rt.genes_nb()); }

#endif
}

} // frictionless
