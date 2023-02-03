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
 * 1. constructor [calling new_transcriptome]
 * 2. new_signature_size
 * 3. init_signature
 * 4. new_swap (on a previously initialized signature)
 * 5. score
 *
 * @note This is only double-checked with asserts, thus in Debug builds.
 */
class FriedmanScore {
    public:
        /** Constructor
         *
         * @param rt Ranked expression table.
         * @param alpha The alpha parameter of the score function.
         */
        FriedmanScore( const Transcriptome& rt, const double alpha=2, const bool print_progress = false);

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
        /** @} */

        /** @name Intialization
         * Values that are going to be updated at each gene swap, but need to be initialized first.
         * @{ */
             bool _has_init_signature;
        /** @} */

        /** @name Updated when one swap a gene:
         * Members part of the incremental evaluation.
         * @{ */

        /** Cache guard. */
        size_t _cached_gene_in;
        /** Cache guard. */
        size_t _cached_gene_out;

        //! Cache data structures involved in swaping genes: R, A, D.
        CacheSwap _swap_cache;

        /** @name Depending on signature size:
         * Remain constant if the number of genes does not change.
         * @{ */

        /** Cache guard. */
        size_t _cached_signature_size;

        /** Squared the number of genes times cubic number of cells.
         * \f[
         *     B_i(G) = 3|G|^2 m_i(m_i+1)^2
         * \f] */
        std::vector<double> B;

        /** Number of genes times squared number of cells.
         * \f[
         *     C_i(G) = |G|m_i(m_i+1)
         * \f] */
        std::vector<double> C;
        /** @} */


        //! Cache data structure that are stable for a given transcritpome: E, F, GG, T, SSR
        CacheTranscriptome _transcriptome_cache;

        // TODO
        // std::vector<double> logpvals;
        // std::vector<double> S_hats;

    #ifndef NDEBUG
        public:
     #else
        protected:
     #endif

        void clear_cache();

        /** Pre-compute constants for a given transcriptome: E, F, GG, T, SSR.
         *
         * To be called again if the transcriptome has been updated.
         * Called by the constructor for now.
         */
        void new_transcriptome(const bool print_progress = false);

        /** Computes squared root of the log chi-squared. */
        double sqrt_logchisq(const double s, const double m) const;

     public:

        /** Pre-compute constants for a given signature size: B, C.
         *
         * To be called again if the signature size changed.
         */
        void new_signature_size(const size_t signature_size);

        /** Returns true if a first signature has been evaluated.
         *
         * Useful to check whether you have to call init_signature.
         */
        bool has_init_signature();

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

        CacheSwap& swap_cache();
        void swap_cache( const CacheSwap& cache);
        void swap_cache( CacheSwap&& cache);
};

} // frictionless
