#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

#include <eo>
#include <mo>

#define CLUTCHLOG_DEFAULT_DEPTH_BUILT_NODEBUG clutchlog::level::progress
#include <clutchlog/clutchlog.h>

#include <frictionless/frictionless.h>
#include <frictionless/transcriptome.h>
#include <frictionless/parser.h>
#include <frictionless/score.h>
#include <frictionless/eval.h>
#include <frictionless/log.h>

#include "frictionlesser.h"

int main(int argc, char* argv[])
{
    // /*
    eoParser argparser(argc, argv);
    // eoState state;

    const std::string exprsfile = argparser.createParam<std::string>("", "exprs",
        "File of the expressions table to be ranked", 'x', "Data").value();

    const std::string ranksfile = argparser.createParam<std::string>("", "ranks",
        "File of the input ranks table", 'r', "Data").value();

    const std::string signatures = argparser.createParam<std::string>("", "signatures",
        "Name of a file containing candidate/starting signatures", 'i', "Data").value();

    // const bool permute = argparser.createParam<bool>(false, "permute",
    //     "Randomly permute the data to get rid of the signal", 'R', "Data").value();

    const size_t genesetsize = argparser.createParam<size_t>(10, "ngenes",
        "Number of genes in the signatures", 'g', "Parameters").value();

    // const double alpha = argparser.createParam<double>(1, "alpha",
    //     "Score adjustment exponent on the number of genes", 'a', "Parameters").value();

    // const double beta = argparser.createParam<double>(2, "beta",
    //     "Exponent on the log of p-values", 'b', "Parameters").value();

    // const bool optimum = argparser.createParam<bool>(true, "optimum",
    //     "Stop search only when having reach a local optimum", 'o', "Stopping Criterion").value();

    unsigned long long seed = argparser.createParam<long>(0, "seed",
        "Seed of the pseudo-random generator (0 = Epoch)", 's', "Misc").value();

    const std::string log_level = argparser.createParam<std::string>("Progress", "log-level",
        "Maximum depth level of logging (Critical<Error<Warning<Progress<Note<Info<Debug<XDebug, default=Progress)", 'l', "Misc").value();

    const std::string log_file = argparser.createParam<std::string>(".*", "log-file",
        "Regexp indicating which source file is allowed logging (default=all)", 'f', "Misc").value();

    const std::string log_func = argparser.createParam<std::string>(".*", "log-func",
        "Regexp indicating which function is allowed logging (default=all)", 'F', "Misc").value();

   const size_t log_depth = argparser.createParam<size_t>(9999, "log-depth",
        "Maximum stack depth above which logging is not allowed (default=no limit)", 'D', "Misc").value();

    const size_t max_errors = argparser.createParam<size_t>(30, "max-errors",
        "Maximum number of errors reported for each check", 'm', "Misc").value();

    const double epsilon = argparser.createParam<double>(1e-10, "epsilon",
        "Precision for floating-point numbers comparison", 'e', "Misc").value();

    make_verbose(argparser);
    make_help(argparser);
// */
    /*
    std::string exprsfile="../tests/neftel_small_OK.csv";
    std::string ranksfile="";
    std::string signatures="";
    // size_t seed=0;
    std::string log_level="Debug";
    size_t max_errors=30;
    */

    frictionless::clutchlog_config(); // common config
    auto& log = clutchlog::logger();
    ASSERT(log.levels().contains(log_level));
    log.threshold(log_level);
    log.depth(log_depth);
    log.file(log_file);
    log.func(log_func);

    if(seed == 0) {
        seed = std::time(nullptr); // Epoch
    }


    frictionless::Transcriptome tr(max_errors);

    std::string filename;
    bool have_exprs;
    if(exprsfile != "" and ranksfile == "") {
        CLUTCHLOG(progress, "Loading an expressions file...");
        have_exprs = true;
        filename = exprsfile;
        frictionless::NeftelExprParser tableparser(max_errors);
        if(filename == "-") {
            tr = tableparser(std::cin);
        } else {
            tr = load_file(filename, max_errors, tableparser);
        }

    } else if(ranksfile != "" and exprsfile == "") {
        CLUTCHLOG(progress, "Loading a ranks file...");
        have_exprs = false;
        filename = ranksfile;
        frictionless::CommonRankParser tableparser(max_errors);
        if(filename == "-") {
            tr = tableparser(std::cin);
        } else {
            tr = load_file(filename, max_errors, tableparser);
        }

    } else {
        EXIT_ON_ERROR(Invalid_Argument, "You have to load either a rank file or an expression file.");
    }

    check_consistency(tr);
    CLUTCHLOG(note, "OK");

    if(have_exprs) {
        // Convert the expressions to ranks.
        CLUTCHLOG(progress, "Compute ranks...");
        frictionless::Transcriptome ranked = frictionless::rank(tr, /* print_progress */true, epsilon);
        CLUTCHLOG(note, "OK");

        CLUTCHCODE(xdebug,
            if( ranked.ranks().size() < 90 and ranked.genes_nb() < 125) {
                CLUTCHLOG(xdebug, "Ranked table:" << std::endl << ranked.as_art(true) );
            }
        );

        CLUTCHLOG(progress, "Check ranks consistency...");
        try {
            ranked.check_ranks(epsilon);
        } catch(const frictionless::DataSumRanks& e) {
            EXIT_ON_ERROR(DataSumRanks, e.what());
        }
        CLUTCHLOG(note, "OK -- checks passed.");

        CLUTCHLOG(progress, "Output...");
        ranked.as_csv(std::cout, "\t");

        CLUTCHLOG(progress, "Done.");
        exit(static_cast<unsigned char>(Error::No_Error));
    }

    CLUTCHLOG(progress, "Check ranks consistency...");
    try {
        tr.check_ranks(epsilon);
    } catch(const frictionless::DataSumRanks& e) {
        EXIT_ON_ERROR(DataSumRanks, e.what());
    }
    CLUTCHLOG(note, "OK -- checks passed.");

    CLUTCHLOG(progress, "Pre-compute Friedman score cache...");
        frictionless::FriedmanScore frs(tr,2, /*print_progress=*/true);
        frs.new_signature_size(genesetsize);
        ASSERT(tr.genes_nb() >= genesetsize);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Pick a random signature...");
        frictionless::Signature signature(tr.genes_nb());

        std::vector<size_t> indices(tr.genes_nb());
        std::iota(indices.begin(), indices.end(), 0);
        std::mt19937 rg(seed);
        std::shuffle(indices.begin(), indices.end(), rg);
        for(size_t i=0; i<genesetsize; ++i) {
            ASSERT( std::find(tr.genes().begin(), tr.genes().end(), i) != tr.genes().end() );
            signature.select(indices[i]);
            CLUTCHLOGD(xdebug, tr.gene_name(i), 1);
        }
        CLUTCHLOG(debug, "Signature of size: " << genesetsize);
        ASSERT(signature.selected.size() == genesetsize);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Instantiate solver...");
        frictionless::EvalFull feval(frs);
        frictionless::EvalSwap peval(frs);

        // Print stuff on stdlog.
        eoOStreamMonitor logger(std::clog);
            moBestFitnessStat<frictionless::Signature> best;
            logger.add(best);
            moCounterStat<frictionless::Signature> counter;
            logger.add(counter);
        // Print stuff every 10 seconds.
        eoTimedMonitor every(10);
            every.add(logger);
        // Continue search until exhaustion of the neighborhood.
        moTrueContinuator<frictionless::Neighbor> until_end;
        moCheckpoint<frictionless::Neighbor> check(until_end);
            check.add(best);    // Update this stat.
            check.add(counter); // Update this state.
            check.add(every);   // Call this monitor.

        // The actual algorithm.
        frictionless::Neighborhood neighborhood;
        // Hill climber, selecting a random solution among the equal-best ones.
        moRandomBestHC<frictionless::Neighbor>
            search(neighborhood, feval, peval, check);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Solver run...");
        feval(signature);
        CLUTCHLOG(note, "Initial signature: " << signature.str());

        search(signature);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Done.");
}
