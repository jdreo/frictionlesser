#include <sstream>
#include <iostream>

#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"
#include "frictionless/parser.h"

#include <catch2/catch_test_macros.hpp>

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {
    GIVEN( "The Zakiev transcriptome parser" ) {
        frictionless::Transcriptome rt(5);
        frictionless::ZakievParser parser(5);

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

            THEN( "the data loads correctly" ) {
                CHECK_NOTHROW( rt = parser(iss) );
                // Print the ranks table as a color map.
                std::cout << rt.format_ranks(true) << std::endl;
            }
        }
        WHEN( "loading data with inconsistent number of levels" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0 1 2\n"
                "Gene_0\n"
                "0 1 2";
            std::istringstream iss(fake);

            THEN( "an error is raised" ) {
                CHECK_THROWS( rt = parser(iss) );
            }
        }
        WHEN( "loading data with inconsistent number of genes" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0\n"
                "0 1 2\n"
                "Gene_1";
            std::istringstream iss(fake);

            THEN( "an error is raised" ) {
                CHECK_THROWS( rt = parser(iss) );
            }
        }
    }
}
