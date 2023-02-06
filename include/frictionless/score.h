#pragma once

#include <vector>

#include "cache.h"
#include "transcriptome.h"
#include "signature.h"

namespace frictionless {

/** Data structure for computing the score, along with some cached data.
 *
 * Methods updating the internal state are expected to honor the following order.
 * I.E. If you call a method, you should have called before the previous ones.
 * 1. constructor
 * 2. new_transcriptome  OR load_transcriptome_cache
 * 3. new_signature_size OR load_size_cache
 * 4. init_signature
 * 5. new_swap (on a previously initialized signature)
 * 6. score
 *
 * @note This is only double-checked with asserts, thus only in Debug builds.
 */
class FriedmanScore {
    public:
        /** Constructor
         *
         * @param rt Ranked expression table.
         * @param alpha The alpha parameter of the score function.
         */
        FriedmanScore( const Transcriptome& rt, const double alpha=2);

#ifndef NDEBUG
    public: // Needs to be public for testing.
#else
    protected:
#endif
        /** @name Parameters:
         * Chosen at instantiation by the user.
         * @{ */

        /** The ranked transcriptome table. */
        const Transcriptome& _transcriptome;

        /** The user-set parameter of the score. */
        const double alpha;

        /** @} Parameters */

        /** @name Cache management
         * @{ */

        //! New signature cache guard.
        bool _has_init_signature;

        //! Swap cache guard.
        size_t _cached_gene_in;

        //! Swap cache guard.
        size_t _cached_gene_out;

        //! Cache data structures involved in swaping genes: R, A, D.
        CacheSwap _swap_cache;

        /** Cache depending on signature size: B, C.
         *
         * Remain constant if the number of genes does not change.
         */
        CacheSize _size_cache;

        //! cache guard.
        bool _has_size_cache;

        //! Cache data structure that are stable for a given transcritpome: E, F, GG, T, SSR
        CacheTranscriptome _transcriptome_cache;

        //! cache guard.
        bool _has_transcriptome_cache;

     public:
        //! Accessor to get the swap cache.
        CacheSwap& swap_cache();

        //! Accessor to set the swap cache by copy.
        void swap_cache( const CacheSwap& cache);

        //! Accessor to set the swap cache by move.
        void swap_cache( CacheSwap&& cache);

        //! Empty all caches.
        void clear_cache();

        /** Load a transcriptome cache from a stream.
         *
         * @warning STREAM SHOULD BE OPENED IN BINARY MODE.
         */
        void load_transcriptome_cache(std::istream& in);

        /** Load a size cache from a stream.
         *
         * @warning STREAM SHOULD BE OPENED IN BINARY MODE.
         */
        void load_size_cache(std::istream& in);

        /** Save the transcriptome cache to a stream.
         *
         * @warning STREAM SHOULD BE OPENED IN BINARY MODE.
         */
        void save_transcriptome_cache(std::ostream& out);

        /** Save the size cache to a stream.
         *
         * @warning STREAM SHOULD BE OPENED IN BINARY MODE.
         */
        void save_size_cache(std::ostream& out);

        /** Returns true if a transcriptome cache is computed or loaded.
         *
         * Useful to check whether you have to call new_transcriptome.
         */
        bool has_transcriptome_cache();

        /** Returns true if a size cache is computed or loaded.
         *
         * Useful to check whether you have to call new_signature_size.
         */
        bool has_size_cache();

        /** Returns true if a first signature has been evaluated.
         *
         * Useful to check whether you have to call init_signature.
         */
        bool has_init_signature();

        /** @} Cache management */

    #ifndef NDEBUG
        public:
     #else
        protected:
     #endif

        /** Computes squared root of the log chi-squared. */
        double sqrt_logchisq(const double s, const double m) const;

    public:
        /** Pre-compute constants for a given transcriptome: E, F, GG, T, SSR.
         *
         * To be called again if the transcriptome has been updated.
         * Called by the constructor for now.
         */
        void new_transcriptome(const bool print_progress = false);

        /** Pre-compute constants for a given signature size: B, C.
         *
         * To be called again if the signature size changed.
         */
        void new_signature_size(const size_t signature_size);

        /** Compute constants from scratch. */
        void init_signature(const Signature& genes);

        /** Pre-compute constants for a single gene-swap: A, D, R.
         *
         * To be called again if two genes have been swapped.
         *
         * @warning Expects the signature to have been previously initialized.
         *          You must ensure to call init_signature at least once before.
         */
        void new_swap(const size_t gene_in, const size_t gene_out);

        /** Global score. */
        double score(const Signature& genes);
};

} // frictionless
