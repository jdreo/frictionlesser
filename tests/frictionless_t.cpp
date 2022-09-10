#include <sstream>
#include <iostream>

#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"
#include "frictionless/parser.h"
#include "frictionless/score.h"

#include <catch2/catch_all.hpp>

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
        WHEN( "loading data with inconsistent number of geneset" ) {
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
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {1, 2, 3, 4, 5};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "the most simple vector of unsorted values" ) {
        std::vector<double> values = {11, 10, 12, 14, 13};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {2, 1, 3, 5, 4};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "the most simple vector of unsorted values with equalities" ) {
        std::vector<double> values = {10, 10, 12, 14, 14};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {1.5, 1.5, 3, 4.5, 4.5};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "a vector of all equals values" ) {
        std::vector<double> values = {0, 0, 0, 0, 0};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {3, 3, 3, 3, 3};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "a realistic vector" ) {
        std::vector<double> values = {0, 0, 0, 0.3, 0.2, 0.1, 0};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
            THEN( "ranks are correct" ) {
                std::vector<double> expected = {2.5, 2.5, 2.5, 7, 6, 5, 2.5};
                REQUIRE( ranks == expected );
            }
        }
    }

    GIVEN( "an empty vector" ) {
        std::vector<double> values = {};
        WHEN( "ranking it" ) {
            std::vector<double> ranks = frictionless::ranks_of(values, 1e-10);
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
            frictionless::Transcriptome rk = frictionless::rank(exprs, false, 1e-10);

            THEN( "Ranks are consistents" ) {
                REQUIRE(rk.check_ranks(1e-10));
            }
        }
    }
    GIVEN( "Another simple expression table" ) {
        const std::string ssv =
            "GENE    S0  S0  S0  S1  S1  S1\n"
            "G0      1   1   1   1   1   1\n"
            "G1      0   1   2   1   1   1\n"
            "G2      1   1   1   0   1   2\n"
            "G3      0   1   2   0   1   2\n";
        std::istringstream iss(ssv);
        frictionless::NeftelExprParser parser(0);
        frictionless::Transcriptome exprs = parser(iss);

        WHEN( "Ranking the table" ) {
            frictionless::Transcriptome rk = frictionless::rank(exprs, false, 1e-10);

            THEN( "Ranks are consistents" ) {
                REQUIRE(rk.check_ranks(1e-10));
            }
        }
    }
}

SCENARIO( "Friedman cache" ) {
    GIVEN( "A very simple ranked table" ) {
        const std::string ssv =
            "GENE    S0  S0  S0  S1  S1  S1\n"
            "G0      2   2   2   2   2   2\n"
            "G1      1   2   3   2   2   2\n"
            "G2      2   2   2   1   2   3\n"
            "G3      1   2   3   1   2   3\n";
        std::istringstream iss(ssv);
        frictionless::CommonRankParser parser(/*max_errors*/0);
        frictionless::Transcriptome rk = parser(iss);

        WHEN( "Computing transcriptome cache" ) {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);

            THEN( "Transcriptome cache for cells number is consistent" ) {
                for(size_t i : {0,1}) { // 2 samples.
                    REQUIRE(frs.E [i] == 144);
                    REQUIRE(frs.F [i] == 12);
                    REQUIRE(frs.GG[i] == 1.5);
                }
            }
            THEN( "Transcriptome cache for squared ranks sums is consistent" ) {
                //              i  j     cells: 0    1    2
                //              v  v            v    v    v
                REQUIRE(frs.SSR[0][0] == 12); // 2² + 2² + 2²
                REQUIRE(frs.SSR[0][1] == 14); // 1² + 2² + 3²
                REQUIRE(frs.SSR[0][2] == 12);
                REQUIRE(frs.SSR[0][3] == 14);
                REQUIRE(frs.SSR[1][0] == 12);
                REQUIRE(frs.SSR[1][1] == 12);
                REQUIRE(frs.SSR[1][2] == 14);
                REQUIRE(frs.SSR[1][3] == 14);
            }
            THEN( "Transcriptome cache for tie-adjustment factors is consistent" ) {
                //            i  j      ranks: 0    1    2
                //            v  v             v    v    v
                REQUIRE(frs.T[0][0] == 27); // 0³ + 3³ + 0³
                REQUIRE(frs.T[0][1] ==  3); // 1³ + 1³ + 1³
                REQUIRE(frs.T[0][2] == 27);
                REQUIRE(frs.T[0][3] ==  3);
                REQUIRE(frs.T[1][0] == 27);
                REQUIRE(frs.T[1][1] == 27);
                REQUIRE(frs.T[1][2] ==  3);
                REQUIRE(frs.T[1][3] ==  3);
            }
        }
        WHEN( "Computing a two-geneset cache" ) {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frs.new_signature_size(2);

            THEN( "Signature size cache is consistent" ) {
                REQUIRE(frs.B[0] == 576); // 3×2²×3(3+1)²
                REQUIRE(frs.B[1] == 576);
                REQUIRE(frs.C[0] ==  24); //   2×3(3+1)
                REQUIRE(frs.C[1] ==  24);
            }
        }
        WHEN( "Computing a three-geneset cache" ) {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frs.new_signature_size(3);

            THEN( "Signature size cache is consistent" ) {
                REQUIRE(frs.B[0] == 1296); // 3×3²×3(3+1)²
                REQUIRE(frs.B[1] == 1296);
                REQUIRE(frs.C[0] ==   36); //   3×3(3+1)
                REQUIRE(frs.C[1] ==   36);
            }
        }
        WHEN( "Computing two-geneset signature cache from scratch") {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frictionless::Signature geneset(4,false); // 4 geneset.
            // Selects two geneset.
            geneset[0] = true;
            geneset[1] = true;
            frs.new_signature_size(2);
            frs.init_signature(geneset);

            THEN( "Init cache for cellular rank sum is consistent" ) {
                //       cell: c
                //             v
                REQUIRE(frs.R[0] == 3); // 2+1
                REQUIRE(frs.R[1] == 4); // 2+2
                REQUIRE(frs.R[2] == 5); // 2+3
                REQUIRE(frs.R[3] == 4); // 2+2
                REQUIRE(frs.R[4] == 4); // 2+2
                REQUIRE(frs.R[5] == 4); // 2+2
            }
            THEN( "Init cache for sum of squared ranks is consistent" ) {
                // 2 samples.
                REQUIRE(frs.A[0] == 600); // 12(3²+4²+5²)
                REQUIRE(frs.A[1] == 576); // 12(3*4²)
            }
            THEN( "Init cache for average sum of gap to average tied ranks number is consistent" ) {
                // 2 samples.
                REQUIRE(frs.D[0] == 12); // 1/(3-1)*(27+3)-3/(3-1)*2
                REQUIRE(frs.D[1] == 24); // 1/(3-1)*(27+27)-3(3-1)*2
            }
        }
    }
}

SCENARIO( "Friedman score" ) {
    GIVEN( "A very simple ranked table" ) {
        const std::string ssv =
            "GENE    S0  S0  S0  S1  S1  S1\n"
            "G0      2   2   2   2   2   2\n"
            "G1      1   2   3   2   2   2\n"
            "G2      2   2   2   1   2   3\n"
            "G3      1   2   3   1   2   3\n";
        std::istringstream iss(ssv);
        frictionless::CommonRankParser parser(/*max_errors*/0);
        frictionless::Transcriptome rk = parser(iss);
        frictionless::FriedmanScore frs(rk, /*alpha*/2);

        WHEN( "Considering a new signature" ) {
            frictionless::Signature geneset(4,false); // 4 geneset.
            // Selects two geneset.
            geneset[0] = true;
            geneset[1] = true;
            double s;

            THEN( "Scoring raises no error" ) {
                // Two selected genes among 4 available.
                CHECK_NOTHROW(frs.new_signature_size(2));
                CHECK_NOTHROW(frs.init_signature(geneset));
                CHECK_NOTHROW(s = frs.score(geneset));
            }
        }
        WHEN( "Comparing raw and partial evaluation" ) {
            frs.new_signature_size(2);
            double s_A, s_B, s_Bp;

            // First two genes.
            frictionless::Signature A(4,false); // 4 geneset.
            A[0] = true;
            A[1] = true;
            frs.init_signature(A);
            s_A  = frs.score(A);

            // First and last gene.
            frictionless::Signature B(4,false); // 4 geneset.
            B[0] = true;
            B[3] = true;
            frs.init_signature(B);
            s_B = frs.score(B);

            // Remove last gene, select second gene,
            // making the same than signature A.
            size_t jout = 3;
            size_t jin = 1;
            B[jout] = 0;
            B[jin] = 1;
            frs.new_swap(jin, jout);
            s_Bp = frs.score(B);

            THEN( "Scores should be the same" ) {
                REQUIRE( s_A  != Catch::Approx(s_B) );
                REQUIRE( s_B  != Catch::Approx(s_Bp) );
                REQUIRE( s_Bp == Catch::Approx(s_A) );
            }
        }
    }
}
