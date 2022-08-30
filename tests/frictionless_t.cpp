#include <sstream>
#include <iostream>

#include "frictionless/frictionless.hpp"

#include <catch2/catch_test_macros.hpp>

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {
    GIVEN( "A ranked transcriptome instance" ) {
        WHEN( "loading consistent data") {
            const std::string fake =
                "Sample_0 Sample_1 Sample_0 Sample_1\n"
                "Gene_0 0   1   1 0\n"
                "# Comment\n"
                "Gene_1 0.5 0.5   0.5 0.5\n"
                "" // Empty line.
                "Gene_2 0 0 1 1\n"
                "Gene_3 0.5 1 0.5 0   \n";
            std::istringstream iss(fake);
            frictionless::RankedTranscriptome* t;

            THEN( "the data loads correctly" ) {
                CHECK_NOTHROW( t = new frictionless::RankedTranscriptome(iss) );
                // Print the ranks table as a color map.
                std::cout << t->ranks_table(true) << std::endl;
            }
        }
        WHEN( "loading data with inconsistent number of levels" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0 1 2\n"
                "Gene_0\n"
                "0 1 2";
            std::istringstream iss(fake);
            frictionless::RankedTranscriptome* t;

            THEN( "an error is raised" ) {
                CHECK_THROWS( t = new frictionless::RankedTranscriptome(iss) );
            }
        }
        WHEN( "loading data with inconsistent number of genes" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0\n"
                "0 1 2\n"
                "Gene_1";
            std::istringstream iss(fake);
            frictionless::RankedTranscriptome* t;

            THEN( "an error is raised" ) {
                CHECK_THROWS( t = new frictionless::RankedTranscriptome(iss) );
            }
        }
    }
}
