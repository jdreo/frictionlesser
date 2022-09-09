#include <numeric>
#include <vector>

#include "frictionless/frictionless.h"

namespace frictionless {

double sum( const std::vector<double>& t)
{
    return std::accumulate(std::begin(t), std::end(t), 0.0);
}

double sum( const std::vector<size_t>& t)
{
    return std::accumulate(std::begin(t), std::end(t), 0.0);
}

size_t sum( const Signature& genes)
{
    return std::accumulate(std::begin(genes), std::end(genes), 0.0);
}

std::vector<double> ranks_of( const std::vector<double>& exprs, const double epsilon )
{
    CLUTCHCODE(xdebug,
        auto& log = clutchlog::logger();
        auto& out = log.out();
        out << "Ranking: ";
        for(auto e : exprs){ out << e << ", ";}
    );
    const short indexing = 1;

    // Find the permutation of exprs for being sorted.
    std::vector<size_t> at( exprs.size() );
    std::iota( std::begin(at), std::end(at), 0 );
    std::sort( std::begin(at), std::end(at),
        [&exprs](const size_t& i, const size_t& j)->bool {
            return exprs[i] < exprs[j]; } );

    // Compute the ranks.
    std::vector<double> ranks(exprs.size());
    size_t i = 0; // Current index in the permutation.
    while( i < at.size() ) {

        size_t j = i+1; // Running index of the streak of equal ranks.
        std::vector<size_t> streak;
        streak.reserve(exprs.size());
        streak.push_back( i+indexing ); // Rank, not array index.
        while( j < at.size() ) {
            // Epsilon equality of consecutive expressions.
            if( std::abs(exprs[at[i]] - exprs[at[j]]) <= epsilon ) {
                streak.push_back( j + indexing ); // Ranks, not array indices.
                j++;
            } else {
                // The streak ended.
                break;
            }
        } // j
        // The average should work for streak of size 1, too.
        double mean = sum(streak) / (double)streak.size();
        for(size_t k = 0; k < streak.size(); ++k) {
            ranks[at[i+k]] = mean;
        } // k

        if( j >= exprs.size() ) {
            break;
        } else {
            i = j;
        }
    } // i
    CLUTCHCODE(xdebug,
        auto& log = clutchlog::logger();
        auto& out = log.out();
        out << " = ";
        for(auto e : ranks){ out << e << ", ";}
        out << std::endl;
    );
    return ranks;
}

} // frictionless
