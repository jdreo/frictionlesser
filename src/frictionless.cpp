#include <cassert> 
#include <fstream>
#include <sstream>
#include <cmath>

#include <exceptions/exceptions.h>

#include "frictionless/frictionless.hpp"

namespace frictionless {

RankedTranscriptome::RankedTranscriptome( std::istream& input )
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

const std::string& RankedTranscriptome::gene(const size_t j) const
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

const size_t& RankedTranscriptome::index_of(const std::string sample_name) const
{
#ifndef NDEBUG
    return _samples.at(sample_name);
#else
    return _samples[sample_name];
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

const std::map<std::string,size_t>& RankedTranscriptome::samples() const
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
    std::string sample_name;
    getline(input, line);
    std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

    // Cell affiliations header.
    size_t icell = 0;
    size_t isample = 0;
    while(ss >> sample_name) {
        if(not _samples.contains(sample_name)) {
            _samples[sample_name] = isample;
            isample++;
        }
        _affiliations.push_back(_samples[sample_name]);
        _cells_in[_samples[sample_name]].push_back(icell);
        icell++;
    }

    // Rows: Gene, then ranked expressions….
#ifndef NDEBUG
    size_t column = 0;
#endif
    std::string gene_name;
    size_t igene = 0;
    while(true) {
        std::string line;
        double rank;
        getline(input, line);
        std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);
        if(!input) {
            // Mainly catch EOF.
            break;

        } else if(line[0] == '#' or line.empty()) {
            // Catch empty lines or commented lines.
            continue;

        } else if(not isdigit(line[0])) {
            // First column is gene name.
            ss >> gene_name;
            _gene_names.push_back(gene_name);
            _genes.push_back(igene);
            igene++;
#ifndef NDEBUG
            // This should be the first column.
            ASSERT(column == 0);
            column = 0;
#endif
            std::vector<double> ranks_row;
            ranks_row.reserve(_affiliations.size());

            while(ss >> rank) {
                ranks_row.push_back(rank);
#ifndef NDEBUG
                column++;
#endif
            }
            _ranks.push_back(ranks_row);

            if(ranks_row.size() != _affiliations.size()
#ifndef NDEBUG
               or column != _affiliations.size()
#endif
            ) {
                CLUTCHLOG(debug,"       _ranks #= " << _ranks.size());
                CLUTCHLOG(debug,"        igene  = " << igene);
                CLUTCHLOG(debug,"       _genes #= " << _genes.size());
                CLUTCHLOG(debug,"        icell  = " << icell);
                CLUTCHLOG(debug,"_affiliations #= " << _affiliations.size());
                CLUTCHLOG(debug,"      isample  = " << isample);
                CLUTCHLOG(debug,"     _samples #= " << _samples.size());
                CLUTCHLOG(debug,"    _cells_in #= " << _cells_in.size());
                std::ostringstream msg;
                msg << "Inconsistent data"
                    << " (ranks_row #=" << ranks_row.size()
                    << " , affiliations #=" << _affiliations.size()
#ifndef NDEBUG
                    << " , column=" << column
#endif
                    << ")";
                RAISE(DataInconsistent, msg.str());
            }
#ifndef NDEBUG
            // End of row => new column.
            column = 0;
#endif
        } else {
            RAISE(DataRowFormat, "Line is neither EOF, starting with #, empty, or not starting with a digit.");
        }
    }

    CLUTCHLOG(progress, "Check ranks table consistency...");
    CLUTCHLOG(debug,"       _ranks #= " << _ranks.size());
    CLUTCHLOG(debug,"        igene  = " << igene);
    CLUTCHLOG(debug,"       _genes #= " << _genes.size());
    CLUTCHLOG(debug,"        icell  = " << icell);
    CLUTCHLOG(debug,"_affiliations #= " << _affiliations.size());
    CLUTCHLOG(debug,"      isample  = " << isample);
    CLUTCHLOG(debug,"     _samples #= " << _samples.size());
    CLUTCHLOG(debug,"    _cells_in #= " << _cells_in.size());
    if(_genes.size() != _ranks.size()
       or      icell != _affiliations.size()
       or    isample != _cells_in.size()
       or    isample != _samples.size()
       or        _ranks.size() == 0
       or _affiliations.size() == 0
       or      _samples.size() == 0
    ) {
        std::ostringstream msg;
        msg << "Inconsistent data "
            << "(ranks shape: " << _ranks.size()
            << " gene rows, " << _genes.size() << " genes"
            << ", " << icell << " cells"
            << ", " << _samples.size() << " samples"
            << ")";
        RAISE(DataInconsistent, msg.str());
    }
    std::ostringstream msg_genes;
    for(size_t j=0; j < _genes.size(); ++j) {
        if(_ranks.at(j).size() == 0 or _ranks.at(j).size() != _genes.size()) {
            msg_genes << "\tGene `" << _gene_names.at(j) << "` has " << _ranks.at(j).size() << " ranks, should be " << _genes.size() << std::endl;
        }
    }
    if(msg_genes.str() != "") {
        RAISE(DataInconsistent, "\n"+msg_genes.str());
    }

    CLUTCHLOG(debug, "Check sum of ranks across cells...");
    std::ostringstream msg_rank_sum;
    for([[maybe_unused]] const auto& [s,i] : this->samples()) {
        for(size_t j : this->genes()) {
            double sum_ranks = 0;
            for(size_t c : this->cells(i)) {
                sum_ranks += this->rank(c,j);
            }
            const double true_sum_ranks = this->cells_nb(i) * (this->cells_nb(i)-1) / 2;
            if(sum_ranks != true_sum_ranks) {
                msg_rank_sum << "\t - Gene `" << this->gene(j)
                    << "` has ranks sum of " << sum_ranks
                    << " but it should be " << true_sum_ranks
                    << std::endl;
            }
        } // for j in genes
    } // for i in samples
    if(msg_rank_sum.str() != "") {
        RAISE(DataSumRanks, "\n"+msg_rank_sum.str());
    }
    CLUTCHLOG(note, "OK");
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
std::string RankedTranscriptome::ranks_table(bool values) const
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
    for([[maybe_unused]] const auto& [s, i] : rt.samples() ) { // All samples.
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
