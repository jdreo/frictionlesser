#pragma once

#include <string>
#include <vector>

#include <eo>
#include <ga/eoBit.h>

#include <exceptions/exceptions.h>
#include <clutchlog/clutchlog.h>

// Make asserts (de)clutchable.
#define ASSERT(EXPR) { CLUTCHFUNC(critical, assert, EXPR); }

namespace frictionless {

EXCEPTION(Exception, DataError);
    EXCEPTION(DataError, DataInconsistent);
    EXCEPTION(DataError, DataRowFormat);
    EXCEPTION(DataError, DataSumRanks);

double sum( const std::vector<double>& t);
double sum( const std::vector<size_t>& t);

std::vector<double> ranks_of( const std::vector<double>& exprs, const double epsilon=1e-6 );


/** The score used to define the quality of a signature is a floating point number.
 *
 * @see Signature
 */
using Score = double;

/** A signature is a vector of boolean.
 *
 * Each bit=1 means "this gene is part of the signature,
 * each bit=0, "this gene is NOT part of the signature".
 */
using Signature = eoBit<Score,bool>;
 // We use Paradiseo/eo/eoBit with a Score as fitness,
 // and vector<bool> as data structure.

} // frictionless
