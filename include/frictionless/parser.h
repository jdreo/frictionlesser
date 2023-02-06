#pragma once

#include <iostream>

#include <frictionless/transcriptome.h>

namespace frictionless {

//! (Ridiculously complex) separator management for parsing files with istreams.
struct delimiter_ctype : std::ctype<char> {

    //! Add characters to the underlying table mask.
    static const mask* make_table(const std::string& delims);

    //! Constructor taking a list of separators.
    delimiter_ctype(const std::string& delims, ::size_t refs = 0);
};

/** Load a transcriptome from the given input stream.
 *
 * Expects tabular, tabs-separated data of the shape:
 *
 * |  XXXX   |  Aff_0  | … |  Aff_i | … |  Aff_m |
 * |---------|---------|---|--------|---|--------|
 * | gene_0  |  r_00   | … |  r_i0  | … |  r_n0  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_j  |  r_0j   | … |  r_ij  | … |  r_mj  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_n  |  r_0n   | … |  r_in  | … |  r_mn  |
 *
 * @note The first cell of the first row may be ignored,
 *       @see Transcriptome
 */
class TranscriptomeParser
{
    public:
        /** The main call interface. */
        virtual Transcriptome operator()(std::istream& input);

        /** Constructor.
         *
         * @param more_separators List of additional columns/lines separators for parsing files (on top of C defaults, defaults to ",;" for parsing CSV).
         * @param ignore_header_first If `true`, will not load the first element of the header row.
         * @param cell_to_sample A regular expression that will be removed from the cell affiliation, in order to create the sample name.
         * @param errors_max_print Maximum number of identical errors to log. If there are more errors than this number, will draw a random sample of errors to display.
         */
        TranscriptomeParser(
            const std::string more_separators = ",;",
            const bool ignore_header_first = false,
            const std::string cell_to_sample = "",
            const size_t errors_max_print = 20
        );

        virtual ~TranscriptomeParser() {}

    protected:
        /** Additional separators for parsing data files. */
        const std::string _more_separators;

        /** Whether or not to ignore the first element of the header. */
        const bool _ignore_header_first;

        /** Regex to transform a cell name in a sample name
         *
         * Is passed to `std::regex_replace`, you may use one capture group.
         */
        const std::string _cell_to_sample;

        /** The parsed table. */
        Transcriptome _rt;

    protected:
        /** Parse the header. Extract samples and cells. */
        size_t load_header(const std::string& line);

        /** Parse a row. Avoid comments and empty lines. */
        size_t load_row(const std::string& line, size_t& igene);

        /** Parse a gene. Extract the gene name and values. */
        size_t load_gene(std::istringstream& ss, size_t& igene);
};


/** Load a transcriptome from the given input stream.
 *
 * This converts Cell names to their corresponding sample
 * by removing their suffix (e.g. "BT771-P8-A01" -> "BT771-P08").
 *
 * @note Suffix is detected with the following regex: "-[A-Z][0-9]{2}$".
 *
 * Expects tabular, tabs-separated data of the shape:
 *
 * |  GENE   |  Cell_0 | … | Cell_i | … | Cell_m |
 * |---------|---------|---|--------|---|--------|
 * | gene_0  |  r_00   | … |  r_i0  | … |  r_n0  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_j  |  r_0j   | … |  r_ij  | … |  r_mj  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_n  |  r_0n   | … |  r_in  | … |  r_mn  |
 *
 * @note The first cell of the first row may be anything,
 *       but a warning will be printed if it is not "GENE".
 */
class NeftelExprParser : public TranscriptomeParser
{
    public:
        /** Constructor.
         *
         * @param errors_max_print Maximum number of identical errors to log. If there are more errors than this number, will draw a random sample of errors to display.
         */
        NeftelExprParser(const size_t errors_max_print = 20);

        virtual ~NeftelExprParser(){}
};


/** Load a ranked transcriptome from the given input stream.
 *
 * @note This expects the header to display cells affiliations
 *       (i.e. their sample), not cells IDs.
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
 * @warning The first cell of the first row should not exists.
 */
class ZakievRankParser : public TranscriptomeParser
{
    public:
        /** Constructor.
         *
         * @param errors_max_print Maximum number of identical errors to log. If there are more errors than this number, will draw a random sample of errors to display.
         */
        ZakievRankParser(const size_t errors_max_print = 20);

        virtual ~ZakievRankParser(){}
};


/** Load a transcriptome from the given input stream.
 *
 * @note This expects the header to display cells affiliations
 *       (i.e. their sample), not cells IDs.
 *
 * Expects tabular, tabs-separated data of the shape:
 *
 * |  GENE   |  Aff_0  | … |  Aff_i | … |  Aff_m |
 * |---------|---------|---|--------|---|--------|
 * | gene_0  |  r_00   | … |  r_i0  | … |  r_n0  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_j  |  r_0j   | … |  r_ij  | … |  r_mj  |
 * |    …    |    …    | … |   …    | … |   …    |
 * | gene_n  |  r_0n   | … |  r_in  | … |  r_mn  |
 *
 * @note The first cell of the first row may be anything,
 *       but a warning will be printed if it is not "GENE".
 */
class CommonRankParser : public TranscriptomeParser
{
    public:
        /** Constructor.
         *
         * @param errors_max_print Maximum number of identical errors to log. If there are more errors than this number, will draw a random sample of errors to display.
         */
        CommonRankParser(const size_t errors_max_print = 20);

        virtual ~CommonRankParser(){}
};

} // frictionless
