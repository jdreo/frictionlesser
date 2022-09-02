#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <clutchlog/clutchlog.h>

#include "frictionless/frictionless.h"
#include "frictionless/parser.h"

namespace frictionless {


TranscriptomeParser::TranscriptomeParser(
    const bool ignore_header_first,
    const std::string cell_to_sample,
    const size_t errors_max_print
) :
    _ignore_header_first(ignore_header_first)
{}


Transcriptome TranscriptomeParser::operator()( std::istream& input )
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

    return _rt;
}

size_t TranscriptomeParser::load_header(const std::string& line)
{
    size_t i = -1;
    // Cell affiliations header.
    std::istringstream iss(line);
    std::string name;
    std::map<std::string,size_t> loaded;

    if(_ignore_header_first) {
        iss >> name;
        if(name != "GENE") {
            CLUTCHLOG(warning, "First element of header was not 'GENE' but `" << name
                << "`. I will not load this field as a sample.");
        }
    }
    while(iss >> name) {
        size_t current;
        // Extract the sample name from the cell name.
        std::string sample_name;
        if(_cell_to_sample.size() > 0) {
            std::regex cell_id(_cell_to_sample);
            sample_name = std::regex_replace(name, cell_id, "");
        } else {
            sample_name = name;
        }
        if(not loaded.contains(sample_name)) {
            current = ++i;
            loaded[sample_name] = current;
            _rt._samples.push_back(current);
            _rt._sample_names.push_back(sample_name);
        } else {
            current = loaded[sample_name];
        }
        _rt._affiliations.push_back(current);
        _rt._cells_in[current].push_back(_rt._cells_nb);
        _rt._cells_nb++;
    }
    ASSERT(loaded.size() == i+1);
    ASSERT(_rt._samples.size() == i+1);
    return loaded.size();
}

long TranscriptomeParser::load_row(const std::string& line, size_t& igene)
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

size_t TranscriptomeParser::load_gene(std::istringstream& ss, size_t& igene)
{
    std::vector<double> ranks_row;
    ranks_row.reserve(_rt._affiliations.size());

    // First ncolumn is gene name.
    std::string gene_name;
    ss >> gene_name;
    _rt._gene_names.push_back(gene_name);
    _rt._genes.push_back(igene);
    igene++;

    double rank;
    size_t ncolumn = 0;
    while(ss >> rank) {
        ranks_row.push_back(rank);
        ncolumn++;
    }
    _rt._ranks.push_back(ranks_row);

    return ncolumn;
}

NeftelParser::NeftelParser(const size_t errors_max_print)
    : TranscriptomeParser(true, "\\w*(-[A-Z][0-9]{2})$", errors_max_print)
{}

ZakievParser::ZakievParser(const size_t errors_max_print)
    : TranscriptomeParser(false, "", errors_max_print)
{}

} // frictionless
