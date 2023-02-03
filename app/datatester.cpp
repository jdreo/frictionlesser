#include <string>
#include <iostream>
#include <fstream>

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
        "Filename of the expressions table to be ranked", 'x', "Data").value();

    const std::string ranksfile = argparser.createParam<std::string>("", "ranks",
        "Filename of the input ranks table", 'r', "Data").value();

    const std::string signatures = argparser.createParam<std::string>("", "signatures",
        "Name of a file containing candidate/starting signatures", 'i', "Data").value();

    // const bool permute = argparser.createParam<bool>(false, "permute",
    //     "Randomly permute the data to get rid of the signal", 'R', "Data").value();

    // const size_t genes = argparser.createParam<size_t>(50, "genes",
    //     "Number of genes in the signatures", 'g', "Parameters").value();

    // const double alpha = argparser.createParam<double>(1, "alpha",
    //     "Score adjustment exponent on the number of genes", 'a', "Parameters").value();

    // const double beta = argparser.createParam<double>(2, "beta",
    //     "Exponent on the log of p-values", 'b', "Parameters").value();

    // const bool optimum = argparser.createParam<bool>(true, "optimum",
    //     "Stop search only when having reach a local optimum", 'o', "Stopping Criterion").value();

    // const long seed = argparser.createParam<long>(0, "seed",
    //     "Seed of the pseudo-random generator (0 = Epoch)", 's', "Misc").value();

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

    frictionless::Transcriptome tr(max_errors);

    std::string filename;
    bool have_exprs;
    if(exprsfile != "" and ranksfile == "") {
        have_exprs = true;
        filename = exprsfile;
        frictionless::NeftelExprParser tableparser(max_errors);
        if(filename == "-") {
            tr = tableparser(std::cin);
        } else {
            tr = load_file(filename, max_errors, tableparser);
        }

    } else if(ranksfile != "" and exprsfile == "") {
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

    CLUTCHLOG(debug, "Test signatures data structures...");
        frictionless::Signature geneset(tr.genes_nb());
        geneset.select(0);
        geneset.select(1);
        CLUTCHCODE(xdebug,
            geneset.printOn(std::clog);
            std::clog << std::endl;

            for(size_t j : geneset.selected) {
                std::clog << tr.gene_name(j) << " ";
            }
            std::clog << std::endl;
        );
    CLUTCHLOG(debug, "OK");

    CLUTCHLOG(progress, "Pre-compute Friedman score cache...");
        frictionless::FriedmanScore frs(tr,2,true);
        const size_t geneset_nb = geneset.selected.size();
        CLUTCHLOG(debug, "Signature of size " << geneset_nb);
        frs.new_signature_size(geneset_nb);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Compute Friedman score from scratch...");
        frs.init_signature(geneset);
        geneset.fitness(frs.score(geneset));
        CLUTCHLOG(debug, geneset);
    CLUTCHLOG(note, "OK");

    for(size_t i=2; i < 4; ++i) {
        CLUTCHLOG(progress, "Swap two genes and update...");
            geneset.reject(i-1);
            geneset.select(i);

            frs.new_swap(i, i-1);
            geneset.fitness(frs.score(geneset));
            CLUTCHLOG(debug, geneset);
        CLUTCHLOG(note, "OK");
    }

    CLUTCHLOG(progress, "Full eval...");
        frictionless::EvalFull feval(frs);
        geneset.invalidate();
        feval(geneset);
        CLUTCHLOG(info, "Fully evaled solution:" << geneset);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Partial eval on neighbor...");
        frictionless::EvalSwap peval(frs);
        frictionless::Neighbor neighbor(geneset.selected.size());
        neighbor.set(1,3);
        peval(geneset, neighbor);
        CLUTCHLOG(info, "Partially evaled neighbor: " << neighbor);
    CLUTCHLOG(progress, "Equivalent full eval...");
        neighbor.move(geneset);
        feval(geneset); // Compare with full eval.
        CLUTCHLOG(info, "Fully evaled solution: " << geneset);
        CLUTCHLOG(debug, "Solution: " << geneset << ", neighbor: " << neighbor);
        CLUTCHLOG(debug, "Solution fitness: " << geneset.fitness() << ", neighbor fitness: " << neighbor.fitness());
        ASSERT(geneset.fitness() == neighbor.fitness());
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Done.");
}

