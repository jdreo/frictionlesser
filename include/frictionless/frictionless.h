#pragma once

#include <string>
#include <vector>

#include <eo>
#include <ga/eoBit.h>

#include <exceptions/exceptions.h>

#include "log.h"
#include "signature.h"

// Make asserts (de)clutchable.
#define ASSERT(EXPR) { CLUTCHFUNC(critical, assert, EXPR); }

namespace frictionless {

EXCEPTION(Exception, DataError);
    EXCEPTION(DataError, DataInconsistent);
    EXCEPTION(DataError, DataRowFormat);
    EXCEPTION(DataError, DataSumRanks);

double sum( const std::vector<double>& t);
double sum( const std::vector<size_t>& t);

/** Compute the rank of the given vector.
 *
 * @param exprs The vector to be ranked.
 * @param epsilon Precision at which one test for equality.
 */
std::vector<double> ranks_of( const std::vector<double>& exprs, const double epsilon );

} // frictionless
