#pragma once

#include <vector>
#include <iostream>

namespace frictionless {

struct Serialize
{
    virtual void save(std::ostream& out) const = 0;
    virtual void load(std::istream& in ) = 0;

    template<typename T>
    void save(std::ostream& out, T scalar, const size_t depth = 0) const
    {
        static_assert(std::is_trivial_v<T>);
        CLUTCHLOGD(xdebug, "save scalar: " << scalar, depth);
        out.write(reinterpret_cast<char*>(&scalar), sizeof(T));
    }

    void save(std::ostream& out, const std::vector<double>& vec, const size_t depth = 0) const
    {
        CLUTCHLOGD(xdebug, "save vector of size: " << vec.size(), depth);
        this->save(out, vec.size(), depth+1);
        for(const double& v : vec) {
            this->save(out, v, depth+1);
        }
    }

    void save(std::ostream& out, const std::vector<std::vector<double>>& table, const size_t depth = 0) const
    {
        CLUTCHLOGD(xdebug, "save table of size: " << table.size(), depth);
        this->save(out, table.size(), depth+1);
        for(const std::vector<double>& vec : table) {
            this->save(out, vec, depth+1);
        }
    }


    template<typename T>
    double load_scalar(std::istream& in, const size_t depth = 0) const
    {
        static_assert(std::is_trivial_v<T>);
        T scalar;
        in.read(reinterpret_cast<char*>(&scalar), sizeof(T));
        CLUTCHLOGD(xdebug, "loaded scalar: " << scalar, depth);
        return scalar;
    }

    std::vector<double> load_vector(std::istream& in, const size_t depth = 0) const
    {
        CLUTCHLOGD(xdebug, "load vector", depth);
        auto size = this->load_scalar<std::vector<double>::size_type>(in);

        std::vector<double> vec;
        vec.reserve(size);
        for(size_t i = 0; i < size; ++i) {
            vec.push_back( this->load_scalar<double>(in, depth+1) );
        }
        return vec;
    }

    std::vector<std::vector<double>> load_table(std::istream& in, const size_t depth = 0) const
    {
        CLUTCHLOGD(xdebug, "load table", depth);
        auto size = this->load_scalar<std::vector<std::vector<double>>::size_type>(in);

        std::vector<std::vector<double>> table;
        table.reserve(size);
        for(size_t i = 0; i < size; ++i) {
            table.push_back( this->load_vector(in, depth+1) );
        }
        return table;
    }
};

class CacheTranscriptome : public Serialize
{
    public:
        virtual void save(std::ostream& out) const override
        {
            Serialize::save(out, E);
            Serialize::save(out, F);
            Serialize::save(out, GG);

            Serialize::save(out, T);
            Serialize::save(out, SSR);
        }

        virtual void load(std::istream& in) override
        {
            E   = Serialize::load_vector(in);
            F   = Serialize::load_vector(in);
            GG  = Serialize::load_vector(in);

            T   = Serialize::load_table(in);
            SSR = Serialize::load_table(in);
        }

    public:

        void reserve(const size_t samples_nb)
        {
            E  .reserve(samples_nb);
            F  .reserve(samples_nb);
            GG .reserve(samples_nb);
            T  .reserve(samples_nb);
            SSR.reserve(samples_nb);
        }

        void clear()
        {
            E  .clear();
            F  .clear();
            GG .clear();
            T  .clear();
            SSR.clear();
        }

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


/** Cache data structure involved in swaping two genes.
 *
 *  It's purpose is to easily move/copy it around objects.
 * FIXME put code in lib.
 */
class CacheSwap
{
    public:
        CacheSwap( // copy
            const std::vector<double>& r,
            const std::vector<double>& a,
            const std::vector<double>& d
            ) : _R(r), _A(a), _D(d), _has_cache(true)
        {}

        CacheSwap( // move
            std::vector<double>&& r,
            std::vector<double>&& a,
            std::vector<double>&& d
            ) : _R(std::move(r)), _A(std::move(a)), _D(std::move(d)), _has_cache(true)
        {}

        CacheSwap() : _has_cache(false) {}

        void reserve(const size_t samples_nb, const size_t cells_nb)
        {
            _A.reserve(samples_nb);
            _D.reserve(samples_nb);
            _R.reserve(cells_nb);
        }

        void clear()
        {
            _R.clear();
            _A.clear();
            _D.clear();
        }

        /** Sum of ranks across genes, for each cell.
         *   FIXME: $c$ or $c_i$?
         * \f[
         *      R_c(G) = \sum_{t\in G} r_{c_i,t}
         * \f] */
        std::vector<double>& R() {return _R;}
        /** @} */

        /** Sum of squared ranks.
        * \f[
        *     A_i(G) = 12\sum_{c=1}^{m_i} R_c(G)^2
        * \f] */
        std::vector<double>& A() {return _A;}

        /** Average sum of gap to average tied ranks number.
         * \f[
         *     D_i(G) = \frac{1}{m_i-1}\sum_{v = 1}^{|G|} \left(\left(\sum_{a=1}^{g_{v}} t_{v,a}^3\right) - m_i\right)
         * \f] */
        std::vector<double>& D() {return _D;}

    protected:
        std::vector<double> _R;
        std::vector<double> _A;
        std::vector<double> _D;

        bool _has_cache;
};


}
