#pragma once

#include <string>
#include <vector>
#include <map>

namespace frictionless {

/** The data strutcure holding the table of ranked expressions.
 */
class Transcriptome {
    public:
        /** Shortcut for the underlying table data structure. */
        using Type = std::vector<std::vector<double>>;

        /** Table of ranks, covering each gene/cell pair.
         *
         * Ranks may be half-ranks, so we use floating-point numbers.
         */
        Type _ranks;

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

        /** All cells indices at once. */
        std::vector<size_t> _cells_all;

        /** The sample names. */
        std::vector<std::string> _sample_names;

        /** The samples indices. */
        std::vector<size_t> _samples;

    public:
        /** Constructor.
         *
         * @param errors_max_print Maximum number of identical errors to log. If there are more errors than this number, will draw a random sample of errors to display.
         */
        Transcriptome( const size_t errors_max_print = 20 );

        /** @name Accessors on ranks.
         * @{ */

        /** Returns the current rank table. */
        const Type& ranks() const;

        /** Returns the ranks row (containing all cells) for the j-th gene. */
        const std::vector<double>& ranks(const size_t j) const;

        /** Returns the rank for the c-th cell and the j-th gene. */
        const double& rank(const size_t c, const size_t j) const;
        double& rank(const size_t c, const size_t j);

        /** @} */
        /** @name Accessors on genes
         * @{ */

        /** Returns the genes indices. */
        const std::vector<size_t>& genes() const;

        /** Returns the total number of genes. */
        size_t genes_nb() const;

        /** Returns the gene names list. */
        const std::vector<std::string>& gene_names() const;

        /** Returns the name of the j-th gene. */
        const std::string& gene_name(const size_t j) const;

        /** @} */
        /** @name Accessors on samples.
         * @{ */

        /** Returns the samples map. */
        const std::vector<size_t>& samples() const;

        /** Returns the number of samples. */
        size_t samples_nb() const;

        /** Returns the sample names. */
        const std::vector<std::string>& sample_names() const;

        /** Returns the i-th sample name. */
        const std::string& sample_name(const size_t i) const;

        /** Returns the affiliations list. */
        const std::vector<size_t>& affiliations() const;

        /** Returns the name of the c-th affiliation. */
        const size_t& affiliation(const size_t c) const;

        /** @} */
        /** @name Accessors on cells.
         * @{ */

        /** Returns the total number of cells. */
        size_t cells_nb() const;

        /** Returns the number of cells in the i-th sample. */
        size_t cells_nb(const size_t i) const;

        /** Returns the indices of all cells across any samples. */
        std::vector<size_t> cells() const;

        /** returns the cells in the i-th sample. */
        const std::vector<size_t>& cells(const size_t i) const;

        /** @} */

        /** Returns an ASCII-art string representing the ranks table. */
        std::string as_art(const bool values = true) const;

        /** Returns an *SV string of the ranks table. */
        std::ostream& as_csv(std::ostream& out = std::cout, const std::string sep="\t") const;

        /** Maximum number of errors to print. */
        size_t _errors_max_print;

        /** Check the table shape. */
        bool check_tables() const;

        /** Check the genes consistency. */
        bool check_genes() const;

        /** Check the ranks consistency. */
        bool check_ranks(const double epsilon) const;

}; // Transcriptome

Transcriptome rank(const Transcriptome& tr, const bool print_progress, const double epsilon);

} // frictionless
