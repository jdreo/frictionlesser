#include <cassert> 
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <set>

#include <exceptions/exceptions.h>

#include "frictionless/frictionless.hpp"

namespace frictionless {

RankedTranscriptome::RankedTranscriptome( std::istream& input,
                                         const size_t errors_max_print )
    : _cells_nb(0), _errors_max_print(errors_max_print)
{
    load(input);
}

const std::vector<std::vector<double>>& RankedTranscriptome::ranks() const
{
    return _ranks;
}

const std::vector<double>& RankedTranscriptome::ranks(const size_t j) const
{
#ifndef NDEBUG
    return _ranks.at(j);
#else
    return _ranks[j];
#endif
}

const double& RankedTranscriptome::rank(const size_t c, const size_t j) const
{
// NOTE the inversion of i and j,
//      to follow the technical report convention.
#ifndef NDEBUG
    return _ranks.at(j).at(c);
#else
    return _ranks[j][c];
#endif

}

const std::vector<std::string>& RankedTranscriptome::gene_names() const
{
    return _gene_names;
}

const std::vector<size_t>& RankedTranscriptome::genes() const
{
    return _genes;
}

size_t RankedTranscriptome::genes_nb() const
{
    return _genes.size();
}

const std::string& RankedTranscriptome::gene_name(
    const size_t j) const
{
#ifndef NDEBUG
    return _gene_names.at(j);
#else
    return _gene_names[j];
#endif
}

const std::vector<size_t>& RankedTranscriptome::affiliations() const
{
    return _affiliations;
}

const size_t& RankedTranscriptome::affiliation(const size_t c) const
{
#ifndef NDEBUG
    return _affiliations.at(c);
#else
    return _affiliations[c];
#endif
}

const std::vector<size_t>& RankedTranscriptome::cells(const size_t i) const
{
#ifndef NDEBUG
    return _cells_in.at(i);
#else
    return _cells_in[i];
#endif
}

size_t RankedTranscriptome::cells_nb(const size_t i) const
{
#ifndef NDEBUG
    return _cells_in.at(i).size();
#else
    return _cells_in[i].size();
#endif
}

size_t RankedTranscriptome::cells_nb() const
{
    return _cells_nb;
}

const std::vector<std::string>& RankedTranscriptome::sample_names() const
{
    return _sample_names;
}

const std::string& RankedTranscriptome::sample_name(const size_t i) const
{
#ifndef NDEBUG
    return _sample_names.at(i);
#else
    return _sample_names[i];
#endif
}

const std::vector<size_t>& RankedTranscriptome::samples() const
{
    return _samples;
}

size_t RankedTranscriptome::samples_nb() const
{
    return _samples.size();
}

void RankedTranscriptome::load( std::istream& input )
{
    std::string line;
    std::getline(input, line);

    size_t nsample = load_header(line);

    CLUTCHLOG(progress, "Parse ranks table data...");
    // Rows: Gene, then ranked expressions….
    size_t igene = 0;
    long ncolumn = 0;
    while(true) {
        if(not input) { // Mainly catch EOF.
            break;
        } else {
            std::string line;
            getline(input, line);
            ncolumn = load_row(line, igene);
        }
    }

    CLUTCHLOG(progress, "Check data consistency...");
    check_tables();
    check_genes();
    check_ranks();
    CLUTCHLOG(note, "OK -- all checks passed.");
}

size_t RankedTranscriptome::load_header(const std::string& line)
{
    size_t i = -1;
    // Cell affiliations header.
    std::istringstream iss(line);
    std::string sample_name;
    std::map<std::string,size_t> loaded;
    while(iss >> sample_name) {
        size_t current;
        if(not loaded.contains(sample_name)) {
            current = ++i;
            loaded[sample_name] = current;
            _samples.push_back(current);
            _sample_names.push_back(sample_name);
        } else {
            current = loaded[sample_name];
        }
        _affiliations.push_back(current);
        _cells_in[current].push_back(_cells_nb);
        _cells_nb++;
    }
    ASSERT(loaded.size() == i+1);
    ASSERT(_samples.size() == i+1);
    return loaded.size();
}

long RankedTranscriptome::load_row(const std::string& line, size_t& igene)
{
    long ncolumn = 0;
    std::istringstream ss(line);

    if(line[0] == '#' or line.empty()) {
        // Catch empty lines or commented lines.
        return ncolumn;

    } else if(not isdigit(line[0])) {
        return load_gene(ss, igene);

    } else {
        RAISE(DataRowFormat, "Line " << igene+1 << " is neither EOF, starting with #, empty, or not starting with a digit.");
    }
}

size_t RankedTranscriptome::load_gene(std::istringstream& ss, size_t& igene)
{
    std::vector<double> ranks_row;
    ranks_row.reserve(_affiliations.size());

    // First ncolumn is gene name.
    std::string gene_name;
    ss >> gene_name;
    _gene_names.push_back(gene_name);
    _genes.push_back(igene);
    igene++;

    double rank;
    size_t ncolumn = 0;
    while(ss >> rank) {
        ranks_row.push_back(rank);
        ncolumn++;
    }
    _ranks.push_back(ranks_row);

    return ncolumn;
}

void RankedTranscriptome::check_tables()
{
    CLUTCHLOG(progress, "Check ranks table consistency...");
    CLUTCHLOG(debug,"       _ranks #= " << _ranks.size());
    // CLUTCHLOG(debug,"        igene  = " << igene);
    CLUTCHLOG(debug,"       _genes #= " << _genes.size());
    CLUTCHLOG(debug,"    _cells_nb  = " << _cells_nb);
    CLUTCHLOG(debug,"_affiliations #= " << _affiliations.size());
    // CLUTCHLOG(debug,"      nsample  = " << nsample);
    CLUTCHLOG(debug,"     _samples #= " << _samples.size());
    CLUTCHLOG(debug,"    _cells_in #= " << _cells_in.size());
    if(     _genes.size() != _ranks.size()
       or       _cells_nb != _affiliations.size()
       or _samples.size() != _cells_in.size()
       or        _ranks.size() == 0
       or _affiliations.size() == 0
       or      _samples.size() == 0
    ) {
        std::ostringstream msg;
        msg << "Inconsistent data "
            << "(ranks shape: " << _ranks.size()
            << " gene rows, " << _genes.size() << " genes"
            << ", " << _cells_nb << " cells"
            << ", " << _samples.size() << " samples"
            << ")";
        RAISE(DataInconsistent, msg.str());
    }
}

void RankedTranscriptome::check_genes()
{
    CLUTCHLOG(progress, "Check all genes consistency...");
    std::vector<std::string> msg_genes;
    for(size_t j=0; j < _genes.size(); ++j) {
        if(   _ranks.at(j).size() == 0
           or _ranks.at(j).size() != _cells_nb
           or _cells_nb != _affiliations.size()
        ) {
            std::ostringstream msg_gene;
            msg_gene << "\tGene `" << _gene_names.at(j)
                     << "`\thas " << _ranks.at(j).size()
                     << " ranks,\tshould be " << _cells_nb
                     << " with " << _affiliations.size() << " affiliations";
            msg_genes.push_back(msg_gene.str());
        }
    }
    if(msg_genes.size() > 0) {
        std::ostringstream msg;
        const size_t max_size = std::min(msg_genes.size(), _errors_max_print);
        std::copy(std::begin(msg_genes), std::begin(msg_genes)+max_size, std::ostream_iterator<std::string>(msg, "\n"));
        if(max_size == _errors_max_print) {
            msg << "[ " << msg_genes.size()-_errors_max_print << " more errors removed ]";}
        RAISE(DataInconsistent, "\n"+msg.str());
    }
}

void RankedTranscriptome::check_ranks()
{
    CLUTCHLOG(debug, "Check sum of ranks across cells...");
    std::vector<std::string> msg_rank_sum;
    for(size_t j : _genes) {
        for(size_t i : _samples) {
            double sum_ranks = 0;
            for(size_t c : this->cells(i)) {
                sum_ranks += this->rank(c,j);
            }
            const double true_sum_ranks = this->cells_nb(i) * (this->cells_nb(i)-1) / 2;
            if(sum_ranks != true_sum_ranks) {
                std::ostringstream msg;
                msg << "\t - Gene " << j << " `" << this->gene_name(j)
                    << "`\thas ranks sum of " << sum_ranks
                    << "\tfor sample " << i
                    << "\tbut it should be " << true_sum_ranks;
               msg_rank_sum.push_back(msg.str());
            }
        } // for j in genes
    } // for i in samples
    if(msg_rank_sum.size() > 0) {
        std::ostringstream msg;
        const size_t max_size = std::min(msg_rank_sum.size(), _errors_max_print);

        // Select some random errors.
        std::vector<std::string> some_msg;
        std::sample(std::begin(msg_rank_sum), std::end(msg_rank_sum),
            std::back_inserter(some_msg),
            max_size, std::mt19937{std::random_device{}()} );

        // Print some messages.
        std::copy(std::begin(some_msg), std::begin(some_msg)+max_size,
            std::ostream_iterator<std::string>(msg, "\n"));
        // Indicate how many were left out.
        if(max_size == _errors_max_print) {
            msg << "\t[ " << msg_rank_sum.size()-_errors_max_print << " more errors not reported ]";}

        RAISE(DataSumRanks, "Sample of " << _errors_max_print << " errors:\n"+msg.str());
    }
}

/** Print a 2D colormap as a 256-color ASCII-art.
 *
 * Convenience function to print the 2D table of transcriptome.
 *
 * Useful to rapidely check or debug the core data structure.
 *
 * @note The y-axis lower cells values are printed as vertical numbers.
 *       The x-axis legend show both min and max of cells.
 *
 * @param values If true, will display one 3-digits value within the colored pixels. Numbers with more than 3 digits are rendered as "+++".
  */
std::string RankedTranscriptome::format_ranks(bool values) const
{
    size_t fg_shift = 128;

    size_t height = _ranks.size();
    ASSERT(height>0);
#ifndef NDEBUG
    size_t width = _ranks.at(0).size();
#else
    size_t width = _ranks[0].size();
#endif
    ASSERT(width>0);

    // Min/max of values.
    auto vmin = static_cast<double>(std::numeric_limits<double>::max());
    auto vmax = static_cast<double>(std::numeric_limits<double>::min());
    for (std::size_t i = 0; i < 0 + height; ++i) {
        for (std::size_t j = 0; j < 0 + width; ++j) {
            const auto val = static_cast<double>(_ranks[i][j]);
            if (not std::isnan(val) and not std::isinf(val)) {
                vmin = std::min(val, vmin);
                vmax = std::max(val, vmax);
    }}}

    const double nb_colors = 255;
    std::ostringstream out;

    // Pixel map
    std::string nan = "╳╳╳";
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            const auto val = static_cast<double>(_ranks[i][j]);
            if (std::isnan(val) or std::isinf(val)) {
                out << nan;
            }
            else{
                auto color = static_cast<size_t>(
                    std::floor((val - vmin) / (vmax - vmin) * (nb_colors - 1)));

                if (color >= nb_colors) {
                    color = static_cast<size_t>(nb_colors - 1);
                }

                if (vmin <= val and val <= vmax) {
                    // Value as background, white or black as foreground.
                    out << "\033[48;2;" << color << ";" << color << ";" << color;
                    size_t fg = (color - fg_shift) % static_cast<size_t>(nb_colors);
                    out << ";38;2;" << fg << ";" << fg << ";" << fg;
                    out << "m";

                    // if asked.
                    if (values /*and j % 2 == 0*/) {
                        std::ostringstream fmt;
                        fmt << val;
                        if (fmt.str().size() == 1) {
                            out << "  " << fmt.str();
                        } else if (fmt.str().size() == 2) {
                            out << " " << fmt.str();
                        } else if (fmt.str().size() == 3) {
                            out << "" << fmt.str();
                        } else {
                            out << "+++";
                        }
                    } else {
                        out << "   ";
                    }

                    out << "\033[0m"; // Colors reset.
                } else {
                    out << nan;
                }
            } // if nan or inf
        } // j
        out << " " << _gene_names[i] << std::endl;
    } // i

    return out.str();
}


FriedmanScore::FriedmanScore( const RankedTranscriptome& rt, const double a)
    : _transcriptome(rt), alpha(a)
{
    ASSERT(rt.ranks( ).size() > 0);
    ASSERT(rt.ranks(0).size() > 0);
    ASSERT(rt.genes_nb() > 0);
    ASSERT(rt.samples_nb() > 0);
    ASSERT(rt.cells_nb(0) > 0);
    ASSERT(alpha > 0);

    const size_t samples_nb = rt.samples_nb();
    A     .reserve(samples_nb);
    D     .reserve(samples_nb);
    B     .reserve(samples_nb);
    C     .reserve(samples_nb);
    E     .reserve(samples_nb);
    F     .reserve(samples_nb);
    GG    .reserve(samples_nb);
    T_ij  .reserve(samples_nb);
    SSR_ij.reserve(samples_nb);

    // Precompute constants.
    for(size_t i : rt.samples() ) { // All samples.
        // Constants depending only on the number of cells
        // in each sample.
        const double m_i = rt.cells_nb(i);
        E .push_back( 3 * m_i * std::pow(m_i+1,2) );
        F .push_back( m_i * (m_i+1) );
        GG.push_back( m_i / (m_i - 1) );

        // Constants for each cells of each sample.
        // FIXME TODO tests
        std::vector<double> SSR_j;
        std::vector<double> T_j;
        for(size_t j : rt.genes()) { // All genes.

            double sum_ranks_sq = 0; // for SSR_ij.
            // Rank => ties count^3, for T_ij.
            // This should work even if ranks are unordered.
            std::map<double,std::vector<double>> ties_cube;

            for(size_t c : rt.cells(i) ) { // All cells of this sample.
                sum_ranks_sq += std::pow(rt.rank(c,j),2);
#ifndef NDEBUG
                auto& t = ties_cube[ rt.rank(c,j) ];
#else
                auto& t = ties_cube.at( rt.rank(c,j) );
#endif
                t.push_back( std::pow(
                    1 + t.size(), 3
                ));

            } // for c in cells
            SSR_j.push_back( sum_ranks_sq );

            // Sum cubic number of ties.
            double sum_ties_cube = 0;
            for([[maybe_unused]] const auto& [r,t] : ties_cube) {
                if(t.size() > 1) { // There were ties for this rank.
                    sum_ties_cube += sum(t);
            }}
#ifndef NDEBUG
            T_j.push_back( sum_ties_cube );
#else
            T_j.push_back( sum_ties_cube );
#endif
        } // for j in genes
        SSR_ij.push_back(SSR_j);
        T_ij.push_back(T_j);
    } // for i in samples

#ifndef NDEBUG
    // Basic size checks.
    ASSERT(E .size() == samples_nb);
    ASSERT(F .size() == samples_nb);
    ASSERT(GG.size() == samples_nb);
    ASSERT(SSR_ij.size() == samples_nb);
    ASSERT(  T_ij.size() == samples_nb);
    for(auto row : SSR_ij) {
        ASSERT(row.size() == rt.genes_nb()); }
    for(auto row : T_ij) {
        ASSERT(row.size() == rt.genes_nb()); }

#endif
}

} // frictionless
