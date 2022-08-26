#pragma once

#include <string>
#include <vector>

#include <eo>
#include <ga/eoBit.h>

#include <clutchlog/clutchlog.h>
// Make asserts (de)clutchable.
// #define ASSERT(LEVEL, ...) { CLUTCHFUNC(LEVEL, assert, __VA_ARGS__) }
#define ASSERT(EXPR, LEVEL) { CLUTCHFUNC(LEVEL, assert, EXPR) }

namespace frictionless {

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
        std::vector<std::string> _genes;

        /** Cell affiliations = index of sample from which this cell was taken. */
        std::vector<size_t> _affiliations;

        /** A map linking index of sample to indices of cells. */
        std::map<size_t,std::vector<size_t>> _cells_in;

        /** A map linking sample name to its index. */
        std::map<std::string,size_t> _samples;

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
        const std::vector<std::vector<double>>& ranks();

        /** Returns the rank row for the i-th gene. */
        const std::vector<double>& ranks(const size_t i);

        /** Returns the genes list. */
        const std::vector<std::string>& genes();

        /** Returns the name of the i-th gene. */
        const std::string& gene(const size_t i);

        /** Returns the affiliations list. */
        const std::vector<size_t>& affiliations();

        /** Returns the name of the i-th affiliation. */
        const size_t& affiliation(const size_t i);

        /** Returns the index of the given sample name. */
        const size_t& index_of(const std::string sample_name);

        /** Returns the number of cells in the given sample index. */
        size_t cells_nb(const size_t sample_id) const;

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


struct Friedman {
    const double alpha;
    // Updated when gene swap.
    std::vector<double> A;
    std::vector<double> D;

    // Constants for a given transcriptome.
    const std::vector<double> B;
    const std::vector<double> C;
    const std::vector<double> E;
    const std::vector<double> F;

    const std::vector<double> GG;
    const std::vector<std::vector<double> > T_ij;
    const std::vector<std::vector<double> > SSR;

    const std::vector<double> Rc;
    const std::vector<double> logpvals;
    const std::vector<double> S_hats;
};

} // frictionless
