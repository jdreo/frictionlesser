#include <cassert> 
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <set>

#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"

namespace frictionless {

Transcriptome::Transcriptome( const size_t errors_max_print )
    : _cells_nb(0), _errors_max_print(errors_max_print)
{
}

const std::vector<std::vector<double>>& Transcriptome::ranks() const
{
    return _ranks;
}

const std::vector<double>& Transcriptome::ranks(const size_t j) const
{
    return _ranks.at(j);
}

const double& Transcriptome::rank(const size_t c, const size_t j) const
{
    return _ranks.at(j).at(c);
}

double& Transcriptome::rank(const size_t c, const size_t j)
{
// NOTE the inversion of i and j,
//      to follow the technical report convention.
#ifndef NDEBUG
    return _ranks.at(j).at(c);
#else
    return _ranks[j][c];
#endif
}

const std::vector<std::string>& Transcriptome::gene_names() const
{
    return _gene_names;
}

const std::vector<size_t>& Transcriptome::genes() const
{
    return _genes;
}

size_t Transcriptome::genes_nb() const
{
    return _genes.size();
}

const std::string& Transcriptome::gene_name(
    const size_t j) const
{
#ifndef NDEBUG
    return _gene_names.at(j);
#else
    return _gene_names[j];
#endif
}

const std::vector<size_t>& Transcriptome::affiliations() const
{
    return _affiliations;
}

const size_t& Transcriptome::affiliation(const size_t c) const
{
    return _affiliations.at(c);
}

const std::vector<size_t>& Transcriptome::cells(const size_t i) const
{
    return _cells_in.at(i);
}

size_t Transcriptome::cells_nb(const size_t i) const
{
    return _cells_in.at(i).size();
}

size_t Transcriptome::cells_nb() const
{
    return _cells_nb;
}

std::vector<size_t> Transcriptome::cells() const
{
    return _cells_all;
}

const std::vector<std::string>& Transcriptome::sample_names() const
{
    return _sample_names;
}

const std::string& Transcriptome::sample_name(const size_t i) const
{
    return _sample_names.at(i);
}

const std::vector<size_t>& Transcriptome::samples() const
{
    return _samples;
}

size_t Transcriptome::samples_nb() const
{
    return _samples.size();
}


bool Transcriptome::check_tables() const
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
    return true;
}

bool Transcriptome::check_genes() const
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
        std::copy(std::begin(msg_genes), std::begin(msg_genes)+max_size,
                  std::ostream_iterator<std::string>(msg, "\n"));
        if(max_size == _errors_max_print) {
            msg << "[ " << msg_genes.size()-_errors_max_print << " more errors removed ]";}
        RAISE(DataInconsistent, "\n"+msg.str());
    }
    return true;
}

bool Transcriptome::check_ranks(const double epsilon) const
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
            if( std::abs(sum_ranks - true_sum_ranks) > epsilon) {
                std::ostringstream msg;
                msg << "\t - Gene " << j << " `" << this->gene_name(j) << "`"
                    << "\tin sample " << i
                    << "\thas ranks sum of " << sum_ranks
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
            msg << "\t[ " << msg_rank_sum.size()-_errors_max_print << " more errors not reported ]";
            RAISE(DataSumRanks, "Sample of " << _errors_max_print << " errors:\n"+msg.str());
        } else {
            RAISE(DataSumRanks, msg_rank_sum.size() << " errors:\n"+msg.str());
        }
    }
    return true;
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
std::string Transcriptome::as_art(bool values) const
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

std::ostream& Transcriptome::as_csv(std::ostream& out, const std::string sep) const
{
    // HEADER
    out << "GENE"; // Following Neftel’s format.
    for(size_t i : samples()) {
        for(size_t c : cells(i)) {
            out << sep << sample_name(i);
        }
    }
    out << "\n";

    // TABLE
    for(size_t j : genes()) {
        out << gene_name(j);
        for(size_t i : samples()) {
            for(size_t c : cells(i)) {
                out << sep << rank(c,j);
            }
        }
        out << "\n";
    }
    return out;
}

Transcriptome rank(const Transcriptome& tr, const bool print_progress)
{
    frictionless::Transcriptome ranked = tr;
    double progress = 0;
    for(size_t j : ranked.genes()) {
        //CLUTCHLOG(xdebug, "Gene "  << j);
        for(size_t i : ranked.samples()) {
            if(print_progress) {
                std::clog << "\r" << static_cast<size_t>(progress / (ranked.genes_nb() * ranked.samples_nb()) * 100) << "%     ";
                std::clog.flush();
            }
            const std::vector<size_t>& cells = ranked.cells(i);
            std::vector<double> exprs;
            exprs.reserve( cells.size() );
            for(size_t c : cells) {
                exprs.push_back( tr.rank(c,j) );
            }
            ASSERT(exprs.size() == cells.size());
            std::vector<double> cell_ranks = frictionless::ranks_of(exprs);
            ASSERT(cell_ranks.size() == exprs.size());
            for(size_t c=0; c < cells.size(); ++c) {
                #ifndef NDEBUG
                    ranked.rank(cells.at(c),j) = cell_ranks.at(c);
                #else
                    ranked.rank(cells[c],j) = cell_ranks[c];
                #endif
            }
            progress++;
        } // for i in samples
    } // for j in genes
    if(print_progress) {
        std::clog << std::endl;
    }
    return ranked;
}

} // frictionless

