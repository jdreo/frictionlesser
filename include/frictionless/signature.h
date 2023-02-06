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

/** The score used to define the fitness (quality) of a signature
 *  is a floating point number along with a cache.
 *
 * The cache is moved along with the fitness, so as to allow for fast partial evaluation.
 *
 * @see Signature
 */
template<class C>
class Score
{
    public:
        //! Scalar type of the score.
        using Type = double;

        //! Attached cache type.
        using Cache = C;

    protected:
        //! The actual score.
        Type _score;

        //! The actual cache.
        Cache _cache;

        //! Cache guard.
        bool _has_cache;

    public:
        //! Empty constructor.
        Score() : _score(std::numeric_limits<Type>::signaling_NaN()), _has_cache(false) {}

        //! Constructor without cache.
        Score(const Type& score) : _score(score), _has_cache(false) {}

        //! Constructor with cache.
        Score(const Type& score, const Cache& cache) :
            _score(score), _cache(cache), _has_cache(true)
        {}

        //! Move constructor.
        Score(Type&& score, Cache&& cache) :
            _score(std::move(score)), _cache(std::move(cache)), _has_cache(true)
        {}

        //! Operator used for silent casting to double.
        operator double() const {return _score;}

        //! Accessor to the actual score.
        Type score() const {return _score;}

        //! Accessor to the cache guard.
        bool has_cache() const {return _has_cache;}

        //! Get a copy of the cache.
        Cache cache() const {ASSERT(_has_cache); return _cache;}

        //! Set the cache by copy.
        void cache(const Cache& c) {_cache = c; _has_cache = true;}

        //! Set the cache by move.
        void cache(Cache&& c) {_cache = std::move(c); _has_cache = true;}

        //! Empty the cache.
        void clear() {
            _cache.clear();
            _has_cache = false;
        }

        //! Serialize the score (not the cache).
        friend std::ostream& operator<<(std::ostream& out, const Score<C>& s)
        {
            // out << s._has_cache;
            // out << " ";
            out << s._score;
            return out;
        }

        //! Deserialize the score (not the cache).
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

//! Proxy name for hiding the default cache type (works only for gene swap).
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
