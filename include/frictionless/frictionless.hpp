#pragma once

#include <string>
#include <vector>

#include <eo>
#include <ga/eoBit.h>

#include <exceptions/exceptions.h>
#include <clutchlog/clutchlog.h>
// Make asserts (de)clutchable.
// #define ASSERT(LEVEL, ...) { CLUTCHFUNC(LEVEL, assert, __VA_ARGS__) }
#define ASSERT(EXPR) { CLUTCHFUNC(critical, assert, EXPR); }

namespace frictionless {

EXCEPTION(Exception, DataError);
    EXCEPTION(DataError, DataInconsistent);
    EXCEPTION(DataError, DataRowFormat);
    EXCEPTION(DataError, DataSumRanks);

template<class T>
double sum( const T& t)
{
    return std::accumulate(std::begin(t), std::end(t), 0);
}


/** The data strutcure holding the table of ranked expressions.
 */
class RankedTranscriptome {
    protected:
        /** Table of ranks, covering each gene/cell pair.
         *
         * Ranks may be half-ranks, so we use floating-point numbers.
         */
        std::vector<std::vector<double>> _ranks;

        /** Gene names. */
        std::vector<std::string> _gene_names;

        /** Genes indices */
        std::vector<size_t> _genes;

        /** Cell affiliations = index of sample from which this cell was taken. */
        std::vector<size_t> _affiliations;

        /** A map linking index of sample to indices of cells. */
        std::map<size_t,std::vector<size_t>> _cells_in;

        /** total number of cells, across all samples. */
        size_t _cells_nb;

        /** A map linking sample name to its index. */
        std::map<std::string,size_t> _samples; // FIXME make it vectors, like for genes.

    public:
        /** Load from a stream of data.
         *
         * Can easily read from either a (ifstream) file
         * or a (istringstream) string.
         *
         * Example with a file:
         * @code
         * std::ifstream ifs;
         * ifs.open("filename");
         * assert(ifs.is_open());
         * frictionless::RankedTranscriptome rt(ifs);
         * ifs.close();
         * @endcode
         *
         * @see load
         */
        RankedTranscriptome(std::istream& input);

        /** Returns the current rank table. */
        const std::vector<std::vector<double>>& ranks() const;

        /** Returns the ranks row (containing all cells) for the j-th gene. */
        const std::vector<double>& ranks(const size_t j) const;

        /** Returns the rank for the c-th cell and the j-th gene. */
        const double& rank(const size_t c, const size_t j) const;

        /** Returns the gene names list. */
        const std::vector<std::string>& gene_names() const;

        /** Returns the genes indices. */
        const std::vector<size_t>& genes() const;

        /** Returns the number of genes. */
        size_t genes_nb() const;

        /** Returns the name of the j-th gene. */
        const std::string& gene(const size_t j) const;

        /** Returns the affiliations list. */
        const std::vector<size_t>& affiliations() const;

        /** Returns the name of the c-th affiliation. */
        const size_t& affiliation(const size_t c) const;

        /** Returns the index of the given sample name. */
        const size_t& index_of(const std::string sample_name) const;

        /** Returns the total number of cells. */
        size_t cells_nb() const;

        /** Returns the number of cells in the i-th sample. */
        size_t cells_nb(const size_t i) const;

        /** returns the cells in the i-th sample. */
        const std::vector<size_t>& cells(const size_t i) const;

        /** Returns the samples map. */
        const std::map<std::string,size_t>& samples() const;

        /** Returns the number of samples. */
        size_t samples_nb() const;

        /** Load a ranked transcriptome from the given input stream.
         *
         * Expects tabular, space-separated data of the shape:
         *
         * ||           Aff_0  | … |  Aff_i | … |  Aff_m |
         * |---------|---------|---|--------|---|--------|
         * | gene_0  |  r_00   | … |  r_i0  | … |  r_n0  |
         * |    …    |    …    | … |   …    | … |   …    |
         * | gene_j  |  r_0j   | … |  r_ij  | … |  r_mj  |
         * |    …    |    …    | … |   …    | … |   …    |
         * | gene_n  |  r_0n   | … |  r_in  | … |  r_mn  |
         *
         * @note The first cell of the first row should not exists.
         */
        void load(std::istream& input);

        /** Returns an ASCII-art string representing the ranks table. */
        std::string ranks_table(const bool values = true) const;

    private:
        size_t load_header(const std::string& line);
        size_t load_gene(std::istringstream& ss, size_t& igene);
        long load_row(const std::string& line, size_t& igene);
        void check_ranks();

}; // RankedTranscriptome

/** The score used to define the quality of a signature is a floating point number.
 *
 * @see Signature
 */
using Score = double;

/** A signature is a vector of boolean.
 *
 * Each bit=1 means "this gene is part of the signature,
 * each bit=0, "this gene is NOT part of the signature".
 */
using Signature = eoBit<Score,bool>;
 // We use Paradiseo/eo/eoBit with a Score as fitness,
 // and vector<bool> as data structure.


class FriedmanScore {
    public:
        FriedmanScore( const RankedTranscriptome& rt, const double alpha=2);

    protected:
        /** @name Parameters:
         * Chosen at instantiation by the user.
         * @{ */

        /** The ranked transcriptome table. */
        const RankedTranscriptome& _transcriptome;

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
