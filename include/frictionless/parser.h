#pragma once

#include <iostream>

#include <frictionless/transcriptome.h>

namespace frictionless {

/** Interface for loading a ranked transcriptome from the given input stream.
 */
class TranscriptomeParser
{
    public:
        virtual Transcriptome operator()(std::istream& input);

        TranscriptomeParser(
            const bool ignore_header_first = false,
            const std::string cell_to_sample = "",
            const size_t errors_max_print = 20
        );

    protected:
        const bool _ignore_header_first;
        const std::string _cell_to_sample;
        Transcriptome _rt;

        size_t load_header(const std::string& line);
        size_t load_gene(std::istringstream& ss, size_t& igene);
        long load_row(const std::string& line, size_t& igene);
};

class NeftelParser : public TranscriptomeParser
{
    public:
        NeftelParser(const size_t errors_max_print = 20);
};

class ZakievParser : public TranscriptomeParser
{
    public:
        ZakievParser(const size_t errors_max_print = 20);
};

} // frictionless
