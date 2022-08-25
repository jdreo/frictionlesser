#pragma once

#include <string>
#include <vector>

namespace frictionless {

class RankedTranscriptome {
    protected:
        /** Gene names. */
        std::vector<std::string> _genes;

        /** Cell affiliations = sample from which this cell was taken. */
        std::vector<std::string> _affiliations;

        /** Table of ranks, covering each gene/cell pair.
         *
         * Ranks may be half-ranks, so we use floating-point numbers.
         */
        std::vector<std::vector<double>> _ranks;

    public:
        /** Load from a stream of data.
         *
         * Can easily read from either a (ifstream) file
         * or a (istringstream) string.
         *
         * Example with a file:
         * @code
         * std::ifstream ifs;
         * ifs.open("filename", std::fstream::in);
         * assert(ifs.is_open());
         * frictionless::RankedTranscriptome rt(ifs);
         * ifs.close();
         * @encode
         *
         * @see load
         */
        RankedTranscriptome(std::istream& input);

        /** Returns the current rank table. */
        const std::vector<std::vector<double>>& ranks();

        /** Load a ranked transcriptome from the given input stream.
         *
         * Expects tabular, space-separated data of the shape:
         *
         * |  Aff_0  | … |  Aff_i | … |  Aff_m |
         * |---------|---|--------|---|--------|
         * | gene_0  |||||
         * |  r_00   | … |  r_i0  | … |  r_n0  |
         * |   …     | … |   …    | … |   …    |
         * | gene_j  |||||
         * |  r_0j   | … |  r_ij  | … |  r_mj  |
         * |   …     | … |   …    | … |   …    |
         * | gene_n  |||||
         * |  r_0n   | … |  r_in  | … |  r_mn  |
         *
         * @note This clears data before loading new ones.
         */
        void load(std::istream& input);

        /** Returns an ASCII-art string representing the ranks table. */
        std::string ranks_table(bool values = true);
}; // RankedTranscriptome


} // frictionless
