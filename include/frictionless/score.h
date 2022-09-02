#pragma once

#include <vector>

#include "frictionless/transcriptome.h"

namespace frictionless {

/** Data structure for computing the score, along with some cached data. */
class FriedmanScore {
    public:
        FriedmanScore( const Transcriptome& rt, const double alpha=2);

    protected:
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

        /** Inverse number of cells (for each sample).
         * \f[
         *    GG_i = \frac{m_i}{m_i-1}
         * \f] */
        std::vector<double> GG;

        /** Tie-adjustment factors for each gene j in each sample i.
         * \f[
         *    T_{ij}=\left(\sum_{a=1}^{g_j} t_{j,a}^3\right)_{i}
         * \f] */
        std::vector<std::vector<double> > T_ij;

        /** Sum of squared ranks of individual gene j in each sample i.
         * \f[
         *     (SSR)_{ij}= \sum_{c=1}^{m_i} r_{cj}^2
         * \f] */
        std::vector<std::vector<double> > SSR_ij;
        /** @} */


        std::vector<double> Rc;
        std::vector<double> logpvals;
        std::vector<double> S_hats;
};

} // frictionless
