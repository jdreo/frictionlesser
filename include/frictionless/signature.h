#pragma once

#include <numeric>
#include <vector>

#include <eo>

#include "frictionless.h"
#include "cache.h"

// To be moved in ParadisEO later on.
#include <frictionless/moBinaryPartition.h>
#include <frictionless/moBinaryPartitionSwapNeighbor.h>
#include <frictionless/moBinaryPartitionSwapNeighborhood.h>

namespace frictionless {

// /** The score used to define the quality of a signature is a floating point number.
//  *
//  * @see Signature
//  */
// using Score = double;

template<class C>
class Score
{
    public:
        using Type = double;
        using Cache = C;

    protected:
        Type _score;
        Cache _cache;
        bool _has_cache;

    public:
        Score() : _score(std::numeric_limits<Type>::signaling_NaN()), _has_cache(false) {}

        Score(const Type& score) : _score(score), _has_cache(false) {}

        Score(const Type& score, const Cache& cache) :
            _score(score), _cache(cache), _has_cache(true)
        {}

        Score(Type&& score, Cache&& cache) :
            _score(std::move(score)), _cache(std::move(cache)), _has_cache(true)
        {}

        operator double() const {return _score;}
        Type score() const {return _score;}

        bool has_cache() const {return _has_cache;}

        Cache cache() const {ASSERT(_has_cache); return _cache;}
        void cache(const Cache& c) {_cache = c; _has_cache = true;}
        void cache(Cache&& c) {_cache = std::move(c); _has_cache = true;}

        void clear() {
            _cache.clear();
            _has_cache = false;
        }

        friend std::ostream& operator<<(std::ostream& out, const Score<C>& s)
        {
            // out << s._has_cache;
            // out << " ";
            out << s._score;
            return out;
        }

        friend std::istream& operator>>(std::istream& in, Score<C>& score)
        {
            Type value;
            // bool has_cache;
            // in >> has_cache;
            in >> value;
            score._score = value;
            score._has_cache = false; // Cache is not serialized.
            // score._has_cache = has_cache;
            return in;
        }
};

using Fitness = Score<CacheSwap>;

/** A signature is a partition of a binary space.
 *
 * Genes' indices are either `selected` or `rejected`.
 */
using Signature = moBinaryPartition<Fitness>;
// We use Paradiseo/mo/moBinaryPartition with Score as a fitness type,
// which use two std::set<size_t> as main data structures.

/** Neighbor of a binary partition that is one swap away. */
using Neighbor = moBinaryPartitionSwapNeighbor<Signature>;

/** Neighborhood that is one swap away. */
using Neighborhood = moBinaryPartitionSwapNeighborhood<Signature>;

} // frictionless
