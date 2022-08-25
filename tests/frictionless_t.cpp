#include <sstream>
#include <iostream>

#include "frictionless/frictionless.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace frictionless;

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {

    GIVEN( "A ranked transcriptome instance" ) {
        frictionless::RankedTranscriptome table;

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

            THEN( "the table loads correctly" ) {
                CHECK( not table.has_data() );
                CHECK_NOTHROW( table.load(iss) );
                CHECK( table.has_data() );
                INFO( table.colormap(true) );
            }
        }
        WHEN( "loading data with inconsistent number of levels" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0\n"
                "0 1 2";
            std::istringstream iss(fake);

            THEN( "an error is raised" ) {
                CHECK_THROWS( table.load(iss) );
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
                CHECK_THROWS( table.load(iss) );
            }
        }
    }
}
