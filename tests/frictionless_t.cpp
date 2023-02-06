#include <sstream>
#include <iostream>

#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"
#include "frictionless/parser.h"
#include "frictionless/score.h"

#include <catch2/catch_all.hpp>

SCENARIO( "Ranked transcriptome data can be loaded", "[data]") {

    // Test several separators.
    std::string s = GENERATE(as<std::string>{}, " ", "  ", "\t", ",", ";");

    GIVEN( "The Zakiev transcriptome parser" ) {
        frictionless::Transcriptome rt(5);
        frictionless::ZakievRankParser parser(5);

        WHEN( "loading consistent data") {
            const std::string fake =
                "Sample_0"+s+"Sample_1"+s+"Sample_0"+s+"Sample_1\n"
                "Gene_0"+s+"0"+s+""+s+""+s+"1"+s+""+s+""+s+"1"+s+"0\n"
                "# Comment\n"
                "Gene_1"+s+"0.5"+s+"0.5"+s+""+s+""+s+"0.5"+s+"0.5\n"
                "" // Empty line.
                "Gene_2"+s+"0"+s+"0"+s+"1"+s+"1\n"
                "Gene_3"+s+"0.5"+s+"1"+s+"0.5"+s+"0"+s+""+s+""+s+"\n";
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
                "Gene_0"+s+"1"+s+"2\n"
                "Gene_0\n"
                "0"+s+"1"+s+"2";
            std::istringstream iss(fake);

            THEN( "an error is raised" ) {
                CHECK_THROWS( rt = parser(iss) );
            }
        }
        WHEN( "loading data with inconsistent number of geneset" ) {
            const std::string fake =
                "Sample_0\n"
                "Gene_0\n"
                "0"+s+"1"+s+"2\n"
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
            frs.new_transcriptome(false);

            THEN( "Transcriptome cache for cells number is consistent" ) {
                for(size_t i : {0,1}) { // 2 samples.
                    REQUIRE(frs._transcriptome_cache.E [i] == 144);
                    REQUIRE(frs._transcriptome_cache.F [i] == 12);
                    REQUIRE(frs._transcriptome_cache.GG[i] == 1.5);
                }
            }
            THEN( "Transcriptome cache for squared ranks sums is consistent" ) {
                //                                   i  j     cells: 0    1    2
                //                                   v  v            v    v    v
                REQUIRE(frs._transcriptome_cache.SSR[0][0] == 12); // 2² + 2² + 2²
                REQUIRE(frs._transcriptome_cache.SSR[0][1] == 14); // 1² + 2² + 3²
                REQUIRE(frs._transcriptome_cache.SSR[0][2] == 12);
                REQUIRE(frs._transcriptome_cache.SSR[0][3] == 14);
                REQUIRE(frs._transcriptome_cache.SSR[1][0] == 12);
                REQUIRE(frs._transcriptome_cache.SSR[1][1] == 12);
                REQUIRE(frs._transcriptome_cache.SSR[1][2] == 14);
                REQUIRE(frs._transcriptome_cache.SSR[1][3] == 14);
            }
            THEN( "Transcriptome cache for tie-adjustment factors is consistent" ) {
                //                                 i  j      ranks: 0    1    2
                //                                 v  v             v    v    v
                REQUIRE(frs._transcriptome_cache.T[0][0] == 27); // 0³ + 3³ + 0³
                REQUIRE(frs._transcriptome_cache.T[0][1] ==  3); // 1³ + 1³ + 1³
                REQUIRE(frs._transcriptome_cache.T[0][2] == 27);
                REQUIRE(frs._transcriptome_cache.T[0][3] ==  3);
                REQUIRE(frs._transcriptome_cache.T[1][0] == 27);
                REQUIRE(frs._transcriptome_cache.T[1][1] == 27);
                REQUIRE(frs._transcriptome_cache.T[1][2] ==  3);
                REQUIRE(frs._transcriptome_cache.T[1][3] ==  3);
            }
        }
        WHEN( "Computing a two-geneset cache" ) {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frs.new_transcriptome(false);
            frs.new_signature_size(2);

            THEN( "Signature size cache is consistent" ) {
                REQUIRE(frs._size_cache.B[0] == 576); // 3×2²×3(3+1)²
                REQUIRE(frs._size_cache.B[1] == 576);
                REQUIRE(frs._size_cache.C[0] ==  24); //   2×3(3+1)
                REQUIRE(frs._size_cache.C[1] ==  24);
            }
        }
        WHEN( "Computing a three-geneset cache" ) {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frs.new_transcriptome(false);
            frs.new_signature_size(3);

            THEN( "Signature size cache is consistent" ) {
                REQUIRE(frs._size_cache.B[0] == 1296); // 3×3²×3(3+1)²
                REQUIRE(frs._size_cache.B[1] == 1296);
                REQUIRE(frs._size_cache.C[0] ==   36); //   3×3(3+1)
                REQUIRE(frs._size_cache.C[1] ==   36);
            }
        }
        WHEN( "Computing two-geneset signature cache from scratch") {
            frictionless::FriedmanScore frs(rk, /*alpha*/2);
            frs.new_transcriptome(false);
            frictionless::Signature geneset(4); // 4 genes.
            // Selects two geneset.
            geneset.select(0);
            geneset.select(1);
            frs.new_signature_size(2);
            frs.init_signature(geneset);

            THEN( "Init cache for cellular rank sum is consistent" ) {
                //       cell: c
                //             v
                REQUIRE(frs._swap_cache.R[0] == 3); // 2+1
                REQUIRE(frs._swap_cache.R[1] == 4); // 2+2
                REQUIRE(frs._swap_cache.R[2] == 5); // 2+3
                REQUIRE(frs._swap_cache.R[3] == 4); // 2+2
                REQUIRE(frs._swap_cache.R[4] == 4); // 2+2
                REQUIRE(frs._swap_cache.R[5] == 4); // 2+2
            }
            THEN( "Init cache for sum of squared ranks is consistent" ) {
                // 2 samples.
                REQUIRE(frs._swap_cache.A[0] == 600); // 12(3²+4²+5²)
                REQUIRE(frs._swap_cache.A[1] == 576); // 12(3*4²)
            }
            THEN( "Init cache for average sum of gap to average tied ranks number is consistent" ) {
                // 2 samples.
                REQUIRE(frs._swap_cache.D[0] == 12); // 1/(3-1)*(27+3)-3/(3-1)*2
                REQUIRE(frs._swap_cache.D[1] == 24); // 1/(3-1)*(27+27)-3(3-1)*2
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
        frs.new_transcriptome(false);

        WHEN( "Considering a new signature" ) {
            frictionless::Signature geneset(4); // 4 genes.
            // Selects two geneset.
            geneset.select(0);
            geneset.select(1);
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
            frictionless::Signature A(4); // 4 genes.
            A.select(0);
            A.select(1);
            frs.init_signature(A);
            s_A  = frs.score(A);

            // First and last gene.
            frictionless::Signature B(4); // 4 genes.
            B.select(0);
            B.select(3);
            frs.init_signature(B);
            s_B = frs.score(B);

            // Remove last gene, select second gene,
            // making the same than signature A.
            B.reject(3);
            B.select(1);
            frs.new_swap(/*in*/1, /*out*/3);
            s_Bp = frs.score(B);

            THEN( "Scores should be the same" ) {
                REQUIRE( s_A  != Catch::Approx(s_B) );
                REQUIRE( s_B  != Catch::Approx(s_Bp) );
                REQUIRE( s_Bp == Catch::Approx(s_A) );
            }
        }
    }
}

SCENARIO("Score cache save/reload is consistent on fake data and low-level API") {

    double i = GENERATE(as<double>{}, 1, 1e2, 1e3);
    double j = GENERATE(as<double>{}, 2e4, 2e5, 2e6);
    double k = GENERATE(as<double>{}, 3e7, 3e8, 3e9);

    GIVEN("A very simple ranked table") {
        std::ostringstream ssv;
        ssv << "GENE S0  S0  S0  S1  S1  S1\n"
            << "G0" << " " << j << " " << j << " " << j << " " << j << " " << j << " 2\n"
            << "G1" << " " << i << " " << j << " " << k << " " << j << " " << j << " 2\n"
            << "G2" << " " << j << " " << j << " " << j << " " << i << " " << j << " 3\n"
            << "G3" << " " << i << " " << j << " " << k << " " << i << " " << j << " 3\n";
        std::istringstream iss(ssv.str());
        frictionless::CommonRankParser parser(/*max_errors*/0);
        frictionless::Transcriptome rk = parser(iss);
        frictionless::FriedmanScore frs1(rk, /*alpha*/2);
        frictionless::FriedmanScore frs2(rk, /*alpha*/2);
        frs1.new_transcriptome(false);

        WHEN("Saving and reload transcriptome cache") {
            std::ofstream ofs("tmp.dat", std::ios::binary | std::ios::trunc);
            frs1._transcriptome_cache.save(ofs);
            ofs.close();

            std::ifstream ifs("tmp.dat", std::ios::binary);
            frs2._transcriptome_cache.load(ifs);
            ifs.close();

            THEN("Data is the same") {
                REQUIRE(frs1._transcriptome_cache.E   == frs2._transcriptome_cache.E);
                REQUIRE(frs1._transcriptome_cache.F   == frs2._transcriptome_cache.F);
                REQUIRE(frs1._transcriptome_cache.GG  == frs2._transcriptome_cache.GG);
                REQUIRE(frs1._transcriptome_cache.T   == frs2._transcriptome_cache.T);
                REQUIRE(frs1._transcriptome_cache.SSR == frs2._transcriptome_cache.SSR);
            }
        }

        WHEN("Saving and reload size cache") {
            frs1.new_signature_size(3);

            std::ofstream ofs("tmp.dat", std::ios::binary);
            frs1._size_cache.save(ofs);
            ofs.close();

            std::ifstream ifs("tmp.dat", std::ios::binary);
            frs2._size_cache.load(ifs);
            ifs.close();

            THEN("Data is the same") {
                REQUIRE(frs1._size_cache.signature_size == frs2._size_cache.signature_size);
                REQUIRE(frs1._size_cache.B == frs2._size_cache.B);
                REQUIRE(frs1._size_cache.C == frs2._size_cache.C);
            }
        }
    }
}

SCENARIO("Score cache save/reload is consistent on realistic data and high-level API") {

    GIVEN("A small ranked table") {
        std::ostringstream ssv;
        // Neftel et al. data, truncated to the 50th column and the 50th row.
        ssv << "GENE    MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1MGH101-P1  MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH101-P1   MGH100-P5   MGH100-P5\n"
                << "A1BG    24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  21\n"
                << "A1BG-AS1    22.5    22.5    22.5    22.5    22.5    22.5    22.5    46  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    47  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    45  22.5    1   2\n"
                << "A1CF    21  21  21  21  21  42  21  21  43  21  46  44  21  21  21  21  21  47  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  45  21  21  21  1.5 1.5\n"
                << "A2M 3   12  3   42  37  27  17  45  14  8   32  7   10  20  15  31  9   19  39  36  22  11  46  44  3   30  3   24  26  38  16  25  35  43  18  41  47  3   28  13  21  29  34  23  40  33  6   21\n"
                << "A2M-AS1 23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    47  1.5 1.5\n"
                << "A2ML1   32  11.5    33.5    36  24  11.5    44  11.5    43  11.5    39  11.5    11.5    11.5    33.5    11.5    41  11.5    11.5    26  30  37  35  46  27  11.5    40  38  11.5    11.5    23  30  11.5    11.5    45  11.5    25  11.5    30  11.5    47  11.5    11.5    28  11.5    42  11.5    1.5 1.5\n"
                << "A2MP1   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "A4GALT  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "A4GNT   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AA06    24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AAAS    22  22  22  22  22  22  22  22  22  22  22  22  22  46  22  22  22  22  44  22  22  22  22  22  22  22  22  22  22  22  22  47  22  22  22  22  22  22  22  22  22  22  22  45  22  22  22  21\n"
                << "AACS    24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  12\n"
                << "AACSP1  23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    47  23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    1.5 1.5\n"
                << "AADAC   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AADACL2 24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AADACL3 24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AADACL4 24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AADAT   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  12\n"
                << "AAED1   23  23  47  23  46  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  1.5 1.5\n"
                << "AAGAB   21  47  21  21  21  21  46  21  21  21  44  21  21  42  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  43  21  21  21  21  45  21  21  21  21  21  21  21  21  21  21  1.5 1.5\n"
                << "AAK1    13  34  12  10  46  3   23  18  17  39  26  21  16  19  32  924 22  14  43  28  2   8   11  31  42  41  37  7   29  40  36  45  33  38  15  30  35  4   44  47  27  25  20  5   6   1   1   2\n"
                << "AAMDC   23  23  23  23  23  23  23  23  23  23  23  23  47  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  46  23  23  23  23  23  23  23  23  23  23  23  23  1.5 1.5\n"
                << "AAMP    17.5    17.5    17.5    17.5    17.5    39  35  41  17.5    47  17.5    36  17.5    17.5    17.5    44  17.5    17.5    42  17.5    17.5    17.5    17.5    45  17.5    38  17.5    17.5    17.5    17.5    40  17.5    17.5    17.5    17.5    17.5    37  17.5    17.5    17.5    17.5    17.5    17.5    46  17.5    17.5    43  21\n"
                << "AANAT   23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  47  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  46  23  23  23  23  23  23  23  1.5 1.5\n"
                << "AAR2    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    45  21.5    46  21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    21.5    44  21.5    43  21.5    21.5    21.5    21.5    21.5    21.5    47  21.5    21.5    21.5    21.5    21.5    21.5    21.5    21\n"
                << "AARD    32  14  33  41  4   2   20  26  3   47  23  30  39  15.5    21  10  28  44  7   6   40  19  25  1   8   24  45  31  46  22  12  15.5    18  36  38  27  37  13  34  43  35  42  17  29  9   5   11  21\n"
                << "AARS    22  22  22  22  22  46  22  22  22  22  45  22  44  22  22  22  22  22  22  22  22  22  47  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  21\n"
                << "AARS2   11  26  40  47  6   9   21  30  35  45  13  23  25  34  3   510 46  28  8   19  4   7   32.5    2   12  22  38  1   43  14  39  15  44  37  27  18  29  31  36  24  41  20  32.5    17  42  16  1   2\n"
                << "AARSD1  22.5    22.5    22.5    22.5    22.5    22.5    22.5    46  22.5    22.5    22.5    22.5    47  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    45  22.5    22.5    22.5    22.5    22.5    22.5    21\n"
                << "AASDH   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "AASDHPPT    22.5    47  22.5    22.5    22.5    22.5    22.5    46  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    45  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    1   2\n"
                << "AASS    9   37  42  28  47  6   8   21  34  27  30  33  20  13  10  712 45  11  4   23  3   31  19  35.5    18  38  41  32  29  15  14  39  22  1   24.5    2   43  24.5    46  40  35.5    26  44  17  16  5   2   1\n"
                << "AATF    21  21  21  21  43  45  21  21  21  21  21  21  21  21  21  21  21  21  42  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  21  46  47  21  21  21  44  21  12\n"
                << "AATK    47  23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    1.5 1.5\n"
                << "AATK-AS1    24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "ABAT    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    45  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    46  22.5    22.5    47  22.5    22.5    22.5    21\n"
                << "ABCA1   25  7   31  7   46  28  22  40  7   36  45  39  7   32  7   42  33  7   43  18  29  14  21  27  20  7   7   7   7   16  35  17  34  44  37  15  7   7   38  19  23  7   47  30  26  41  24  21\n"
                << "ABCA10  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  12\n"
                << "ABCA11P 23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  47  23  23  23  23  23  23  23  23  46  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  23  21\n"
                << "ABCA12  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "ABCA13  6   46  22  43  6   6   25  16  6   30.5    33  13.5    45  36  6   641 34.5    37  6   39  40  20.5    6   27  18.5    38  29  47  23.5    13.5    32  15  42  34.5    30.5    6   6   17  6   23.5    12  18.5    44  20.5    28  26  1   2\n"
                << "ABCA17P 24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "ABCA2   23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    47  23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    23.5    1.5 1.5\n"
                << "ABCA3   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  12\n"
                << "ABCA4   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n"
                << "ABCA5   32  42  37.5    13.5    13.5    13.5    13.5    13.5    29  39  13.5    13.5    13.5    13.5    13.5    13.5    40  13.5    13.5    13.5    13.5    13.5    35  35  13.5    13.5    44  37.5    45  13.5    13.5    33  30  13.5    13.5    31  13.5    41  13.5    46  35  13.5    27  43  13.5    28  47  12\n"
                << "ABCA6   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  12\n"
                << "ABCA7   22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    46  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    45  22.5    22.5    22.5    47  22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    22.5    1.5 1.5\n"
                << "ABCA8   24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  1.5 1.5\n";


        std::istringstream iss(ssv.str());
        frictionless::CommonRankParser parser(/*max_errors*/0);
        frictionless::Transcriptome rk = parser(iss);
        frictionless::FriedmanScore frs1(rk, /*alpha*/2);
        frictionless::FriedmanScore frs2(rk, /*alpha*/2);
        frs1.new_transcriptome(false);

        WHEN("Saving and reload transcriptome cache") {
            std::ofstream ofs("tmp_trans.dat", std::ios::binary | std::ios::trunc);
            frs1.save_transcriptome_cache(ofs);
            ofs.close();

            std::ifstream ifs("tmp_trans.dat", std::ios::binary);
            frs2.load_transcriptome_cache(ifs);
            ifs.close();

            THEN("Data is the same") {
                REQUIRE(frs1._transcriptome_cache.E   == frs2._transcriptome_cache.E);
                REQUIRE(frs1._transcriptome_cache.F   == frs2._transcriptome_cache.F);
                REQUIRE(frs1._transcriptome_cache.GG  == frs2._transcriptome_cache.GG);
                REQUIRE(frs1._transcriptome_cache.T   == frs2._transcriptome_cache.T);
                REQUIRE(frs1._transcriptome_cache.SSR == frs2._transcriptome_cache.SSR);
            }
        }

        WHEN("Saving and reload size cache") {
            frs1.new_signature_size(3);

            std::ofstream ofs("tmp_size.dat", std::ios::binary | std::ios::trunc);
            frs1.save_size_cache(ofs);
            ofs.close();

            std::ifstream ifs("tmp_size.dat", std::ios::binary);
            frs2.load_size_cache(ifs);
            ifs.close();

            THEN("Data is the same") {
                REQUIRE(frs1._size_cache.signature_size == frs2._size_cache.signature_size);
                REQUIRE(frs1._size_cache.B == frs2._size_cache.B);
                REQUIRE(frs1._size_cache.C == frs2._size_cache.C);
            }
        }
    }
}
