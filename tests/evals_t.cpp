#include "frictionless/frictionless.h"
#include "frictionless/transcriptome.h"
#include "frictionless/parser.h"
#include "frictionless/score.h"
#include <frictionless/eval.h>

#include <catch2/catch_all.hpp>

SCENARIO( "Evaluators are consistent" ) {
    GIVEN( "A very simple problem" ) {
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
        frictionless::EvalFull feval(frs);
        frictionless::EvalSwap peval(frs);

        WHEN( "Computing a full evaluation" ) {
            frictionless::Signature geneset(rk.genes_nb());
            geneset.select(0);
            geneset.select(1);
            feval(geneset);

            THEN( "Score is consistent" ) {
                REQUIRE( geneset.selected.size() == 2);
                REQUIRE( not geneset.invalid() );
                REQUIRE( geneset.fitness() == Catch::Approx(0.5) );
            }
        }
        WHEN( "Computing partial evaluations of known signatures" ) {
            frictionless::EvalTest teval(feval);

            frictionless::Signature geneset(rk.genes_nb());
            geneset.select(0);
            geneset.select(1);
            feval(geneset); // We need a first eval.

            frictionless::Neighborhood phood;
            frictionless::Neighbor part;
            phood.init(geneset, part);

            frictionless::Neighborhood fhood;
            frictionless::Neighbor full;
            fhood.init(geneset, full);

            // Use the neighborhood to generate signatures to test.
            size_t i = 0;
            do {
                peval(geneset, part);
                teval(geneset, full);

                DYNAMIC_SECTION( "Neighbor `" << part << "` have the same score with full evaluation" ) {
                    REQUIRE( not part.invalid() );
                    REQUIRE( not full.invalid() );
                    REQUIRE( part.fitness() == Catch::Approx(full.fitness()) );
                }

                phood.next(geneset, part);
                fhood.next(geneset, full);
                i++;
            } while(phood.cont(geneset) and fhood.cont(geneset));
        }
    }
}
