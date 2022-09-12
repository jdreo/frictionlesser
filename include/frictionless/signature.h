#pragma once

#include <eo>


#include <frictionless/moBinaryPartition.h>

namespace frictionless {

/** The score used to define the quality of a signature is a floating point number.
 *
 * @see Signature
 */
using Score = double;

/** A signature is a partition of a binary space.
 *
 * Genes' indices are either `selected` or `rejected`.
 */
using Signature = moBinaryPartition<double>;
// We use Paradiseo/mo/moBinaryPartition with Score as a fitness type,
// which use two std::set<size_t> as main data structures.

} // frictionless
