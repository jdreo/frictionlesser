#include <sstream>
#include <iostream>

#include "frictionless/frictionless.hpp"

#include <catch2/catch_test_macros.hpp>

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {
    GIVEN( "A ranked transcriptome instance" ) {
        WHEN( "loading consistent data") {
            const std::string fake =
                "Sample_0 Sample_1 Sample_0\n"
                "Gene_0\n"
                "0   0   0\n"
                "# Comment\n"
                "Gene_1\n"
                "0.5 1   3\n"
                ""
                "Gene_2\n"
                "0.5 2   2\n"
                "Gene_3\n"
                "1   3   1\n";
            std::istringstream iss(fake);
            frictionless::RankedTranscriptome* t;

            THEN( "the data loads correctly" ) {
                CHECK_NOTHROW( t = new frictionless::RankedTranscriptome(iss) );
                // Print the ranks table as a color map,
                // if test fails.
                INFO( t->ranks_table(true) );
            }
        }
        WHEN( "loading data with inconsistent number of levels" ) {
            const std::string fake =
                "Sample_0\n"
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
