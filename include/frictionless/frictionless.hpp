#pragma once

#include <string>
#include <vector>

namespace frictionless {

class RankedTranscriptome {
    protected:
        /** Protection flag. */
        bool _has_data;

        /** Gene names. */
        std::vector<std::string> _genes;

        /** Cell affiliations. */
        std::vector<std::string> _cells;

        /** Table of ranks, covering each gene/cell pair.
         *
         * Ranks may be half-ranks, so we use floating-point numbers.
         */
        std::vector<std::vector<double>> _ranks;

    public:
        /** Load from a filename.
         *
         * @see load
         */
        RankedTranscriptome(const std::string filename);

        /** Instantiate an empty object.
         *
         * Use this if you want to load data later from a string.
         * @see load
         */
        RankedTranscriptome();

        /** Return true if this instance has already loaded consistent data. */
        bool has_data();

        /** Load a ranked transcriptome from the given input stream.
         *
         * Expects tabular, space-separated data of the shape:
         *
         * | Cell_0  | … | Cell_i | … | Cell_m |
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
}; // RankedTranscriptome

} // frictionless
