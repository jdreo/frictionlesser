#pragma once

#include <string>
#include <sstream>
#include <iostream>

#include "frictionless/parser.h"

namespace frictionless {

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
class ZakievRanksParser : public TranscriptomeParser
{
    public:
        ZakievRanksParser(const size_t errors_max_print = 20)
            : _rt(errors_max_print)
        {}

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
         * frictionless::Transcriptome rt(ifs);
         * ifs.close();
         * @endcode
         */
        virtual Transcriptome operator()(std::istream& input) override;

    protected:
        Transcriptome _rt;

        size_t load_header(const std::string& line);
        size_t load_gene(std::istringstream& ss, size_t& igene);
        long load_row(const std::string& line, size_t& igene);
};

} // frictionless
