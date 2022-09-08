#include <sstream>
#include <iostream>

#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"
#include "frictionless/parser.h"

#include <catch2/catch_test_macros.hpp>

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {
    GIVEN( "The Zakiev transcriptome parser" ) {
        frictionless::Transcriptome rt(5);
        frictionless::ZakievRankParser parser(5);

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
                std::cout << rt.as_art(true) << std::endl;
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

SCENARIO( "Ranking a vector" ) {

    GIVEN( "the most simple vector of sorted values" ) {
        std::vector<double> values = {10, 11, 12, 13, 14};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {0, 1, 2, 3, 4};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "the most simple vector of unsorted values" ) {
        std::vector<double> values = {11, 10, 12, 14, 13};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {1, 0, 2, 4, 3};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "the most simple vector of unsorted values with equalities" ) {
        std::vector<double> values = {10, 10, 12, 14, 14};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {0.5, 0.5, 2, 3.5, 3.5};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "a vector of all equals values" ) {
        std::vector<double> values = {0, 0, 0, 0, 0};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {2, 2, 2, 2, 2};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "a realistic vector" ) {
        std::vector<double> values = {0, 0, 0, 0.3, 0.2, 0.1, 0};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {1.5, 1.5, 1.5, 6, 5, 4, 1.5};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "an empty vector" ) {
        std::vector<double> values = {};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values);
            THEN( "we get an empty vector" ) {
                REQUIRE( ranks.size() == 0 );
            }
        }
    }
}

SCENARIO( "Header cell-to-sample regex replacement" ) {
    GIVEN( "a realistic header" ) {
        std::vector<std::string> header = {
            "BT1160-P1-A01",
            "BT749-P10-B02",
            "MGH100-P5-C03",
            "MGH104negP2-D04",
            "MGH105A-E05",
            "MGH105A-P1-F06",
            "MGH106CD3posP1-G07",
            "MGH110-CD45-P1-H08",
            "MGH122-CD45pos-P1-I09",
            "MGH124-CD45negP1-J10",
            "MGH128CD45neg-P3-K11",
            "MGH66-P03-M99"
        };
        WHEN( "substituting cells for samples" ) {

            std::regex c2s("-[A-Z][0-9]{2}$");

            std::vector<std::string> expected = {
                "BT1160-P1",
                "BT749-P10",
                "MGH100-P5",
                "MGH104negP2",
                "MGH105A",
                "MGH105A-P1",
                "MGH106CD3posP1",
                "MGH110-CD45-P1",
                "MGH122-CD45pos-P1",
                "MGH124-CD45negP1",
                "MGH128CD45neg-P3",
                "MGH66-P03"
            };
            std::vector<std::string> results;
            for(std::string name : header) {
                results.push_back( std::regex_replace(name, c2s, "") );
            }
            THEN( "sample names are correct" ) {
                REQUIRE( results == expected );
            }
        }
    }
}

SCENARIO( "Ranking an expressions table" ) {
    GIVEN( "A simple expression table" ) {
        const std::string ssv =
            "GENE S01-A01 S01-B02 S03-A01 S01-B03 S02-A01 S02-A02\n"
            "A    0       0       0       0       0       0\n"
            "B    0.1     0.2     0       0.3     0.1     0.2\n"
            "C    1.1     1.2     0       1.3     1.1     1.2\n"
            "D    0.3     0.2     1       0.1     0       0\n";
        std::istringstream iss(ssv);
        frictionless::NeftelExprParser parser(0);
        frictionless::Transcriptome exprs = parser(iss);

        WHEN( "Ranking the table" ) {
            frictionless::Transcriptome rk = frictionless::rank(exprs, false);

            THEN( "Ranks are consistents" ) {
                REQUIRE(rk.check_ranks());
            }
        }
    }
}
