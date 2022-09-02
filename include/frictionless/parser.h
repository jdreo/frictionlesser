#pragma once

#include <iostream>

#include <frictionless/transcriptome.h>

namespace frictionless {

/** Interface for loading a ranked transcriptome from the given input stream.
 */
class TranscriptomeParser
{
    public:
        virtual Transcriptome operator()(std::istream& input) = 0;
};

} // frictionless
