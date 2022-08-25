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
    return this->_ranks;
}

void RankedTranscriptome::load( std::istream& input )
{
    // Clear everything, in case it has been filled before.
    _ranks.clear();
    _genes.clear();
    _affiliations.clear();

    std::string line;
    std::string charbuf;
    getline(input, line);
    std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

    // Cell affiliations header.
    while(ss >> charbuf) { _affiliations.push_back(charbuf); }

    // Rows: Gene, then ranked expressions….
#ifndef NDEBUG
    size_t column = 0;
#endif
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
            ss >> charbuf;
            _genes.push_back(charbuf);
#ifndef NDEBUG
            // This should be the first column.
            assert(column == 0);
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
                std::ostringstream msg;
                msg << "Inconsistent data"
                    << " (ranks_row.size=" << ranks_row.size()
                    << " , _affiliations.size=" << _affiliations.size()
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
    if(_genes.size() != _ranks.size()) {
        std::ostringstream msg;
        msg << "Inconsistent data"
            << " (_genes.size=" << _genes.size()
            << " , _ranks.size=" << _ranks.size()
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
std::string RankedTranscriptome::ranks_table(bool values)
{
    size_t fg_shift = 128;

    size_t height = _ranks.size();
    assert(height>0);
    size_t width = _ranks.at(0).size();
    assert(width>0);

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
