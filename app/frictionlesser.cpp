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
    eoParser argparser(argc, argv);
    // eoState state; // TODO

    // Add common options.

    // std::string exprsfile;
    // std::string ranksfile;
    // std::string cachetransfile;
    // std::string cachesizefile;
    // bool cacheonly;
    // double alpha;
    // std::string log_level;
    // std::string log_file;
    // std::string log_func;
    // size_t log_depth;
    // size_t max_errors;
    // double epsilon;

    // std::tie(
    //     exprsfile,
    //     ranksfile,
    //     cachetransfile,
    //     cachesizefile,
    //     cacheonly,
    //     alpha,
    //     log_level,
    //     log_file,
    //     log_func,
    //     log_depth,
    //     max_errors,
    //     epsilon
    // ) = frictionless::make_parser(argparser);

    const std::string exprsfile = argparser.createParam<std::string>("", "exprs",
        "File of the expressions table to be ranked", 'x', "Data").value();

    const std::string ranksfile = argparser.createParam<std::string>("", "ranks",
        "File of the input ranks table", 'r', "Data").value();


    const std::string cachetransfile = argparser.createParam<std::string>("", "cache-transcriptome",
        "File storing the transcriptome-dependant cache", 'T', "Cache").value();

    const std::string cachesizefile = argparser.createParam<std::string>("", "cache-size",
        "File storing the geneset-size-dependant cache", 'S', "Cache").value();

    const bool cacheonly = argparser.createParam<bool>(false, "cache-only",
        "Exit after creating transcriptome and size cache files", 'O', "Cache").value();

    // const std::string signatures = argparser.createParam<std::string>("", "signatures",
        // "Name of a file containing candidate/starting signatures", 'i', "Data").value(); // TODO

    // const bool permute = argparser.createParam<bool>(false, "permute",
    //     "Randomly permute the data to get rid of the signal", 'R', "Data").value(); // TODO

    const double alpha = argparser.createParam<double>(0, "alpha",
        "Score adjustment exponent on the number of genes", 'a', "Parameters").value();

    // const double beta = argparser.createParam<double>(2, "beta",
    //     "Exponent on the log of p-values", 'b', "Parameters").value(); // TODO?

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

    const size_t genesetsize = argparser.createParam<size_t>(10, "ngenes",
        "Number of genes in the signatures", 'g', "Parameters").value();

    unsigned long long seed = argparser.createParam<long>(0, "seed",
        "Seed of the pseudo-random generator (0 = Epoch)", 's', "Misc").value();

    const std::string save_sol = argparser.createParam<std::string>("", "save-sol",
        "File to which save every solution encountered (default='', do not save)", 'o', "Misc").value();

    const bool randsign = argparser.createParam<bool>(false, "random-signature",
        "Do not run the search algorithm but stop after evaluating the first random signature", '1', "Parameters").value();

    /**************************************************************************
     * Messages management.
     *************************************************************************/

    make_verbose(argparser);
    make_help(argparser);

    frictionless::clutchlog_config(); // common config
    auto& log = clutchlog::logger();
    ASSERT(log.levels().contains(log_level));
    log.threshold(log_level);
    log.depth(log_depth);
    log.file(log_file);
    log.func(log_func);


    if(cacheonly and cachetransfile == "") {
        EXIT_ON_ERROR(Invalid_Argument, "You asked for computing cache file(s) but did not indicate which transcriptome cache you want, configure at least --cache-transcriptome.");
    }

    if(seed == 0) {
        seed = std::time(nullptr); // Epoch
    }


    /**************************************************************************
     * Input data management.
     *************************************************************************/

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

    frictionless::check_consistency(tr);
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
    ASSERT(tr.genes_nb() >= genesetsize);
    CLUTCHLOG(note, "OK -- checks passed.");


    /**************************************************************************
     * Cache(s) management.
     *************************************************************************/

    frictionless::FriedmanScore frs(tr, alpha);

    // TRANSCRIPTOME CACHE.

    if(cachetransfile != "" and std::filesystem::exists(cachetransfile)) {
        CLUTCHLOG(progress, "Load Friedman score transcriptome cache...");
        CLUTCHLOG(note,"From file: " << cachetransfile);
        // Mandatory binary, or shit will happen.
        std::ifstream ifs(cachetransfile, std::ios::binary);
        if(ifs.fail()) {
            EXIT_ON_ERROR(Unreadable, "Input transcriptome cache file cannot be read."); }
        ASSERT(ifs.is_open());
        frs.load_transcriptome_cache(ifs);
        ifs.close();
    } else {
        CLUTCHLOG(progress, "Compute Friedman score transcriptome cache...");
        frs.new_transcriptome(/*print_progress=*/true);
        if(cachetransfile != "") {
            CLUTCHLOG(note, "Save cache to file: " << cachetransfile);
            // Mandatory binary, or shit will happen.
            std::ofstream ofs(cachetransfile, std::ios::binary | std::ios::trunc);
            if(ofs.fail()) {
                EXIT_ON_ERROR(Unreadable, "Output transcriptome cache file cannot be wrote."); }
            ASSERT(ofs.is_open());
            frs.save_transcriptome_cache(ofs);
        } else {
            CLUTCHLOG(note, "I will not save the transcriptome cache.");
        }
    }
    ASSERT(frs.has_transcriptome_cache());
    CLUTCHLOG(note, "OK");

    if(cacheonly and cachesizefile == "") {
        CLUTCHLOG(progress, "Asked for transcriptome file cache only, stop here.");
        ASSERT(std::filesystem::exists(cachetransfile));
        EXIT_ON_ERROR(No_Error,"Done");
    }


    // SIZE CACHE.

    if(cachesizefile != "" and std::filesystem::exists(cachesizefile)) {
        CLUTCHLOG(progress, "Load Friedman score size cache...");
        CLUTCHLOG(note,"From file: " << cachesizefile);
        // Mandatory binary, or shit will happen.
        std::ifstream ifs(cachesizefile, std::ios::binary);
        if(ifs.fail()) {
            EXIT_ON_ERROR(Unreadable, "Input transcriptome cache file cannot be read."); }
        ASSERT(ifs.is_open());
        frs.load_size_cache(ifs);
        ifs.close();
    } else {
        CLUTCHLOG(progress, "Compute Friedman score size cache...");
        frs.new_signature_size(genesetsize);
        if(cachesizefile != "") {
            CLUTCHLOG(note, "Save size cache to file: " << cachesizefile);
            // Mandatory binary, or shit will happen.
            std::ofstream ofs(cachesizefile, std::ios::binary | std::ios::trunc);
            if(ofs.fail()) {
                EXIT_ON_ERROR(Unreadable, "Output size cache file cannot be wrote."); }
            ASSERT(ofs.is_open());
            frs.save_size_cache(ofs);
        } else {
            CLUTCHLOG(note, "I will not save the size cache.");
        }
    }
    ASSERT(frs.has_size_cache());
    CLUTCHLOG(note, "OK");

    if(cacheonly) {
        CLUTCHLOG(progress, "Asked for file caches only, stop here.");
        ASSERT(std::filesystem::exists(cachetransfile));
        ASSERT(std::filesystem::exists(cachesizefile));
        EXIT_ON_ERROR(No_Error,"Done");
    }


    /**************************************************************************
     * Search Algorithm.
     *************************************************************************/

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

        // Actually compute the init cache.
        frs.init_signature(signature);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Instantiate solver...");
        eo::rng.reseed(seed);
        frictionless::EvalFull feval(frs);
        frictionless::EvalSwap peval(frs);

        // Print stuff on stdlog.
        eoOStreamMonitor logger(std::clog);
            moBestFitnessStat<frictionless::Signature> best;
            logger.add(best);
            // moCounterStat<frictionless::Signature> counter;
            // logger.add(counter);
        // Print stuff every 10 seconds.
        eoTimedMonitor every(10);
            every.add(logger);
        // Continue search until exhaustion of the neighborhood.
        moTrueContinuator<frictionless::Neighbor> until_end;
        moCheckpoint<frictionless::Neighbor> check(until_end);
            // check.add(counter); // Update the counter state.
            check.add(best);    // Update the best state.
            check.add(every);   // Call this monitor.

        // Save every solutions in a file.
        // Declare here for keeping the scope, but will only be added if asked.
        eoFileMonitor save(save_sol,
            /*delimiter*/" ", /*keep_existing*/true, /*header*/true, /*overwrite*/false);
            moSolutionStat<frictionless::Signature> every_sol;
        if(save_sol != "") {
            save.add(every_sol);  // Output the solution to the file.
            check.add(every_sol); // Update the solution state.
            check.add(save);      // Call the saver.
        }

        // The actual algorithm.
        frictionless::Neighborhood neighborhood;
        // Hill climber, selecting a random solution among the equal-best ones.
        moRandomBestHC<frictionless::Neighbor>
            search(neighborhood, feval, peval, check);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Evaluate first random signature...");
    feval(signature);
        CLUTCHLOG(note, "Initial signature: " << signature.str());
    CLUTCHLOG(note, "OK");

    if(not randsign) {
        CLUTCHLOG(progress, "Solver run...");
            search(signature);
        CLUTCHLOG(note, "OK");
    }

    CLUTCHLOG(progress, "Found signature:");
    CLUTCHLOG(note, signature.str() );

    // Real output.
    std::clog << "# score size genes..." << std::endl;
    std::cout << signature.fitness() << " " << signature.selected.size();
    for(auto i : signature.selected) {
        std::cout << " " << tr.gene_name(i);
    }
    std::cout << std::endl;

    CLUTCHLOG(progress, "Done.");
}
