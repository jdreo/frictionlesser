#include <cassert>
#include <fstream>
#include <sstream>

#include "frictionless/frictionless.hpp"

namespace frictionless {

RankedTranscriptome::RankedTranscriptome(const std::string filename)
    : _has_data(false)
{
    std::ifstream ifs;
    ifs.open(filename, std::fstream::in);
    assert(ifs.is_open());
    load(ifs);
    ifs.close();
}

RankedTranscriptome::RankedTranscriptome()
    : _has_data(false)
{
    // Do nothing.
}

bool RankedTranscriptome::has_data()
{
    return this->_has_data;
}

void RankedTranscriptome::load( std::istream& input )
{
    // Clear everything, in case it has been filled before.
    _ranks.clear();
    _genes.clear();
    _cells.clear();
    _has_data = false;

    std::string line;
    std::string charbuf;
    getline(input, line);
    std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

    // Cell affiliations header.
    while(ss >> charbuf) { _cells.push_back(charbuf); }

    // Rows: Gene, then ranked expressions….
#ifndef NDEBUG
    size_t column=0;
#endif
    while(true) {
        std::string line;
        double rank; 
        getline(input, line);
        std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);
        if(!input) {
            // Mainly catch EOF.
            break;

        } else if(line[0] == '#' || line.empty()) {
            // Catch empty lines or commented lines.
            continue;

        } else if(!isdigit(line[0])) {
            // First column is gene name.
            ss >> charbuf;
            _genes.push_back(charbuf);
#ifndef NDEBUG
            // This should be the first column.
            assert(column == 0);
            column = 0;
#endif

        } else {
            std::vector<double> ranks_row;
            ranks_row.reserve(_cells.size());

            while(ss >> rank) {
                ranks_row.push_back(rank);
#ifndef NDEBUG
                column++;
#endif
            }
            _ranks.push_back(ranks_row);

            if(ranks_row.size() != _cells.size()
#ifndef NDEBUG
               or column != _cells.size()
#endif
            ) {
                std::ostringstream msg;
                msg << "Inconsistent data"
                    << " (ranks_row.size=" << ranks_row.size()
                    << " , _cells.size=" << _cells.size()
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

    _has_data = true;
}

}
