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


std::vector<double> ranks_of( const std::vector<double>& exprs )
{
    // Find the permutation of exprs for being sorted.
    std::vector<size_t> at( exprs.size() );
    std::iota( std::begin(at), std::end(at), 0 );
    std::sort( std::begin(at), std::end(at),
        [&exprs](const size_t& i, const size_t& j)->bool {
            return exprs[i] < exprs[j]; } );

    std::vector<double> ranks(exprs.size());
    size_t i = 0;
    while(i < exprs.size()) {

        size_t j = i+1;
        std::vector<size_t> streak;
        streak.reserve(exprs.size());
        streak.push_back(i);
        for(; j < at.size(); ++j) {
            if( exprs[at[j]] == exprs[at[i]] ) {
                streak.push_back(j);
            } else {
                break;
            }
        } // j
        double mean = sum(streak) / (double)streak.size();
        for(size_t k = 0; k < streak.size(); ++k) {
            ranks[at[i+k]] = mean;
        } // k

        if(j >= exprs.size()) {
            break;
        } else {
            i = j;
        }
    } // i
    return ranks;
}

} // frictionless
