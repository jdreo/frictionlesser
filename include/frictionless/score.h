#pragma once

#include <vector>

#include "frictionless/transcriptome.h"

namespace frictionless {

/** Data structure for computing the score, along with some cached data. */
class FriedmanScore {
    public:
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
        /** @} */

        /** @name Updated when one swap a gene:
         * Members part of the incremental evaluation.
         * @{ */

        /** Cache guard. */
        size_t _cached_gene_in;
        /** Cache guard. */
        size_t _cached_gene_out;

        /** Sum of squared ranks.
        * \f[
        *     A_i(G) = 12\sum_{c=1}^{m_i} R_c(G)^2
        * \f] */
        std::vector<double> A;

        /** Average sum of gap to average tied ranks number.
         * \f[
         *     D_i(G) = \frac{1}{m_i-1}\sum_{v = 1}^{|G|} \left(\left(\sum_{a=1}^{g_{v}} t_{v,a}^3\right) - m_i\right)
         * \f] */
        std::vector<double> D;
        /** @} */


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


        /** @name Constants:
         * Remain constant for a given transcriptome.
         * @{ */

        /** Cubic number of cells (for each sample).
         * \f[
         *     E_i = 3 m_i(m_i+1)^2
         * \f] */
        std::vector<double> E;

        /** Squared number of cells (for each sample).
         * \f[
         *    F_i = m_i(m_i+1)
         * \f] */
        std::vector<double> F;

        /** Inverse numb er of cells (for each sample).
         * \f[
         *    GG_i = \frac{m_i}{m_i-1}
         * \f] */
        std::vector<double> GG;

        /** Tie-adjustment factors for each gene j in each sample i.
         * \f[
         *    T_{ij}=\left(\sum_{a=1}^{g_j} t_{j,a}^3\right)_{i}
         * \f] */
        std::vector<std::vector<double> > T;

        /** Sum of squared ranks of individual gene j in each sample i.
         * \f[
         *     (SSR)_{ij}= \sum_{c=1}^{m_i} r_{cj}^2
         * \f] */
        std::vector<std::vector<double> > SSR;
        /** @} */


        /** @name Intialization:
         * Values that are going to be updated at each gene swap, but need to be initialized first.
         * @{ */

        /** Sum of ranks across genes, for each cell.
         * \f[
         *      R_c(G) = \sum_{t\in G} r_{c_i,t}
         * \f] */
        std::vector<double> R;

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
        void new_transcriptome();

        /** Computes squared root of the log chi-squared. */
        double sqrt_logchisq(const double s, const double m) const;

     public:

        /** Pre-compute constants for a given signature size: B, C.
         *
         * To be called again if the signature size changed.
         */
        void new_signature_size(const size_t signature_size);

        /** Compute constants from scratch. */
        void init_signature(Signature genes);

        /** Pre-compute constants for a single gene-swap: A, D.
         *
         * To be called again if two genes have been swapped.
         */
        void new_swap(const size_t gene_in, const size_t gene_out);

        /** Global score. */
        double score(Signature genes);

};

} // frictionless
