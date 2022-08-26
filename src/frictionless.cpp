#include <cassert> 
#include <fstream>
#include <sstream>
#include <cmath>

#include "frictionless/frictionless.hpp"

namespace frictionless {

RankedTranscriptome::RankedTranscriptome( std::istream& input )
{
    load(input);
}

const std::vector<std::vector<double>>& RankedTranscriptome::ranks()
{
    return _ranks;
}

const std::vector<double>& RankedTranscriptome::ranks(const size_t i)
{
#ifndef NDEBUG
    return _ranks.at(i);
#else
    return _ranks[i];
#endif
}

const std::vector<std::string>& RankedTranscriptome::genes()
{
    return _genes;
}

const std::string& RankedTranscriptome::gene(const size_t i)
{
#ifndef NDEBUG
    return _genes.at(i);
#else
    return _genes[i];
#endif
}

const std::vector<size_t>& RankedTranscriptome::affiliations()
{
    return _affiliations;
}

const size_t& RankedTranscriptome::affiliation(const size_t i)
{
#ifndef NDEBUG
    return _affiliations.at(i);
#else
    return _affiliations[i];
#endif
}

const size_t& RankedTranscriptome::index_of(const std::string sample_name)
{
#ifndef NDEBUG
    return _samples.at(sample_name);
#else
    return _samples[sample_name];
#endif
}

size_t RankedTranscriptome::cells_nb(const size_t sample_id) const
{
#ifndef NDEBUG
    return _cells_in.at(sample_id).size();
#else
    return _cells_in[sample_id].size();
#endif
}

size_t RankedTranscriptome::samples_nb() const
{
    return _samples.size();
}

void RankedTranscriptome::load( std::istream& input )
{
    // Clear everything, in case it has been filled before.
    // _ranks.clear();
    // _genes.clear();
    // _affiliations.clear();
    // _cells_in.clear();
    // _samples.clear();

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
            _genes.push_back(gene_name);
#ifndef NDEBUG
            // This should be the first column.
            ASSERT(column == 0, error);
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
                throw std::runtime_error(msg.str());
            }
#ifndef NDEBUG
            // End of row => new column.
            column = 0;
#endif
        } else {
            throw std::runtime_error("Line is neither EOF, starting with #, empty, or not starting with a digit.");
        }
    }

    CLUTCHLOG(debug,"       _ranks #= " << _ranks.size());
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
    ) {
        std::ostringstream msg;
        msg << "Inconsistent data "
            << "(ranks shape: " << _ranks.size() << "×"  << _ranks.at(0).size()
            << ", "  << _genes.size() << " genes"
            << ", " << icell << " cells"
            << ", " << _samples.size() << " samples"
            << ")";
        throw std::runtime_error(msg.str());
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
std::string RankedTranscriptome::ranks_table(bool values) const
{
    size_t fg_shift = 128;

    size_t height = _ranks.size();
    ASSERT(height>0, error);
#ifndef NDEBUG
    size_t width = _ranks.at(0).size();
#else
    size_t width = _ranks[0].size();
#endif
    ASSERT(width>0, error);

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
        out << " " << _genes[i] << std::endl;
    } // i

    return out.str();
}
}
