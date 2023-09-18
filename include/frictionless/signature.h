#pragma once

#include <numeric>
#include <vector>
#include <utility>

#include <eo>

#include "frictionless.h"
#include "cache.h"

// To be moved in ParadisEO later on.
#include <frictionless/moBinaryPartition.h>
#include <frictionless/moBinaryPartitionSwapNeighbor.h>
#include <frictionless/moBinaryPartitionSwapNeighborhood.h>

namespace frictionless {

// FIXME move the code in the src/ dir.
struct ScoreDetails
{
    using ValueType = double;
    using VectorType = std::vector<ValueType>;
    ValueType value;
    VectorType values_by_samples;

    ScoreDetails(const ValueType& v, const VectorType& vs) : value(v), values_by_samples(vs)
    { }

    bool operator==(const ScoreDetails& other)
    {
        if(this->value != other.value) {return false;}
        if(this->values_by_samples.size() != other.values_by_samples.size()) {return false;}
        for(size_t i=0; i < this->values_by_samples.size(); ++i) {
            if(this->values_by_samples[i] != other.values_by_samples[i]) {return false;}
        }
        return true;
    }

    //! Serialize the score.
    friend std::ostream& operator<<(std::ostream& out, const ScoreDetails& s)
    {
        out << s.value;
        out << " " << s.values_by_samples.size();
        for(auto v : s.values_by_samples) {
            out << " " << v;
        }
        return out;
    }

    //! Deserialize the score.
    friend std::istream& operator>>(std::istream& in, ScoreDetails& score)
    {
        double value;
        size_t size;
        std::vector<double> values_by_samples;
        in >> value;

        in >> size;
        values_by_samples.reserve(size);
        for(size_t i=0; i<size; ++i) {
            double v;
            in >> v;
            values_by_samples.push_back(v);
        }
        ASSERT(values_by_samples.size() == size);
        score = ScoreDetails(value, values_by_samples);
        return in;
    }
};

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
        // using Type = double;

        /** Type returned by the "global score" function.
         *
         * Holds both the global score and the Friedman scores by samples.
         */
        using Type = ScoreDetails;

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
        Score() : _score(std::numeric_limits<Type::ValueType>::signaling_NaN(), Type::VectorType()), _has_cache(false) {}

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

        //! Copy constructor.
        // Score( const Score& s ) :
        //     _score(s._score.value, s._score.values_by_samples), _cache(s._cache), _has_cache(s._has_cache)
        // {}

        //! Operator used for silent casting to double.
        operator double() const {return _score.value;}

        //! Accessor to the actual score.
        Type::ValueType score() const {return _score.value;}

        //! Accessor to the actual score.
        Type score_details() const {return _score;}

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
            // bool has_cache;
            // in >> has_cache;
            in >> score._score;

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
