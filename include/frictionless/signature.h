#pragma once

#include <eo>


#include <frictionless/moBinaryPartition.h>
#include <frictionless/moBinaryPartitionSwapNeighbor.h>
#include <frictionless/moBinaryPartitionSwapNeighborhood.h>

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

/** Neighbor of a binary partition that is one swap away. */
using Neighbor = moBinaryPartitionSwapNeighbor<Signature>;

/** Neighborhood that is one swap away. */
using Neighborhood = moBinaryPartitionSwapNeighborhood<Signature>;

} // frictionless
