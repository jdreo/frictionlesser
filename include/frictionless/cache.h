#pragma once

#include <vector>
#include <iostream>

#include "log.h"

namespace frictionless {

/** Basic serialization capabilities.
 *
 * Helper class that provides raw serialization of scalars,
 * vector of doubles and tables of doubles.
 *
 * @warning EXPECT STREAMS IN BINARY MODE.
 *          If saved cache fails to reload, you probably saved/opened in the default text mode.
 *
 * A derived class is expected to implement `save` and `load`,
 * calling the protected methods on the members it wants to cache,
 * while paying attention to the order of their serialization.
 *
 * @note The type of serialization used here is raw memory dump.
 *       While it is fast and compact, it is not suitable for sharing data.
 *       It is only suitable for internal use in the very same executable.
 *       If the cached files are to be read by another binary version,
 *       another software or on another machine, it will probably fail.
 *
 * You are responsible to ensure that the serialization calls
 * follow the same order in `save` and `load`.
 */
class Serialize
{
    public:
        //! Interface for calling a cache write.
        virtual void save(std::ostream& out) const = 0;

        //! Interface for calling a cache read.
        virtual void load(std::istream& in ) = 0;

    protected:
        //! Serialize any scalar.
        template<typename T>
        void save(std::ostream& out, T scalar, const size_t depth = 0) const
        {
            static_assert(std::is_trivial_v<T>);
            CLUTCHLOGD(xdebug, "save scalar: " << scalar, depth);
            out.write(reinterpret_cast<char*>(&scalar), sizeof(T));
        }

        //! Serialize a vector of double.
        void save(std::ostream& out, const std::vector<double>& vec, const size_t depth = 0) const;

        /** Serialize a table of double (i.e. a vector of vector of double).
         *
         * @note This does not check whether each row fas the same length.
         */
        void save(std::ostream& out, const std::vector<std::vector<double>>& table, const size_t depth = 0) const;

        //! Deserialize any scalar.
        template<typename T>
        double load_scalar(std::istream& in, const size_t depth = 0) const
        {
            static_assert(std::is_trivial_v<T>);
            T scalar;
            in.read(reinterpret_cast<char*>(&scalar), sizeof(T));
            CLUTCHLOGD(xdebug, "loaded scalar: " << scalar, depth);
            return scalar;
        }

        //! Deserialize a vector of double.
        std::vector<double> load_vector(std::istream& in, const size_t depth = 0) const;

        /** Deserialize a table of double (i.e. a vector of vector of double).
         *
         * @note This does not check whether each row fas the same length.
         */
        std::vector<std::vector<double>> load_table(std::istream& in, const size_t depth = 0) const;
};

/** Cache for Friedman score's intermediate results that are tied to a given transcriptome.
 *
 * That is:
 *     - vectors: E, F, GG;
 *     - tables: T, SSR
 */
class CacheTranscriptome : public Serialize
{
    public:
        //! Serialize.
        virtual void save(std::ostream& out) const override;

        //! Deserialize.
        virtual void load(std::istream& in) override;

    public:

        //! Reserve memory for all members.
        void reserve(const size_t samples_nb);

        //! Empty all cache.
        void clear();

    public:
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

};

/** Cache for Friedman score's intermediate results that are tied to a given signature size.
 *
 * That is vectors B and C.
 */
class CacheSize : public Serialize
{
    public:
        //! Serialize.
        virtual void save(std::ostream& out) const override;

        //! Deserialize.
        virtual void load(std::istream& in) override;

    public:

        //! Reserve memory for all members.
        void reserve(const size_t samples_nb);

        //! Empty all cache.
        void clear();

        /** @name Depending on signature size:
         * Remain constant if the number of genes does not change.
         * @{ */

        /** Cache guard. */
        size_t signature_size;

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
};


/** Cache data structure involved in swaping two genes.
 *
 *  It's purpose is to easily move/copy it around objects.
 */
class CacheSwap
{
    public:
        //! Copy constructor.
        CacheSwap( // copy
            const std::vector<double>& r,
            const std::vector<double>& a,
            const std::vector<double>& d
        );

        //! Move constructor.
        CacheSwap( // move
            std::vector<double>&& r,
            std::vector<double>&& a,
            std::vector<double>&& d
        );

        CacheSwap();

        //! Reserve memory for all members.
        void reserve(const size_t samples_nb, const size_t cells_nb);

        //! Empty all cache.
        void clear();

        /** Sum of ranks across genes, for each cell.
         *   FIXME: $c$ or $c_i$?
         * \f[
         *      R_c(G) = \sum_{t\in G} r_{c_i,t}
         * \f] */
        std::vector<double> R;

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

    protected:
        //! Cache guard.
        bool _has_cache;
};

}
