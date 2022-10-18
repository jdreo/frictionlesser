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

//! Error codes returned on exit.
enum class Error : unsigned char {
    No_Error = 0,
    No_File = 2, // ENOENT
    Invalid_Argument = 22, // EINVAL
    Unreadable = 77, // EBADFD
    Missing_Argument = 132,
    DataInconsistent,
    DataRowFormat,
    DataSumRanks,
    Unknown = 255
};

//! Macro to hide the necessary verbose casting.
#ifdef WITH_CLUTCHLOG
    #define EXIT_ON_ERROR(err_code, msg) { \
        if(static_cast<unsigned char>(Error::err_code) > 0) { \
            CLUTCHLOG(critical, "CRITICAL ERROR"); \
            CLUTCHCODE(critical, \
                std::cerr << msg << std::endl; \
            ); \
        } \
        exit(static_cast<unsigned char>(Error::err_code)); \
    }
#else
    #define EXIT_ON_ERROR(err_code, msg) { \
        if(static_cast<unsigned char>(Error::err_code) > 0) { \
            std::cerr << "CRITICAL ERROR: " << msg << std::endl; \
        } \
        exit(static_cast<unsigned char>(Error::err_code)); \
    }
#endif


static frictionless::Transcriptome load_file(
    const std::string filename,
    const size_t max_errors,
    frictionless::TranscriptomeParser& tableparser)
{

    if(not std::filesystem::exists(filename)) {
        EXIT_ON_ERROR(No_File, "Input data file does not exists."); }

    CLUTCHLOG(progress, "Open `" << filename << "`...");

    std::ifstream ifs(filename);
    if(ifs.fail()) {
        EXIT_ON_ERROR(Unreadable, "Input data file cannot be read."); }
    ASSERT(ifs.is_open());

    frictionless::Transcriptome tr(max_errors);
    // frictionless::NeftelExprParser tableparser(max_errors);
    tr = tableparser(ifs);
    ifs.close();

    CLUTCHLOG(progress, "Loaded "
        << tr.genes().size() << " genes, "
        << tr.affiliations().size() << " cells, and "
        << tr.samples_nb() << " samples.");

    CLUTCHCODE(xdebug,
        if( tr.ranks().size() < 90 and tr.genes_nb() < 125) {
            CLUTCHLOG(xdebug, "Loaded table:" << std::endl << tr.as_art(true) );
        }
    );

    return tr;
}

static void check_consistency(const frictionless::Transcriptome& tr)
{
    CLUTCHCODE(debug,
        CLUTCHLOG(debug, "Samples:");
        for(const auto& i : tr.samples()) {
            CLUTCHLOG(debug, "    - " << tr.sample_name(i) << ": " << tr.cells(i).size() << " cells" );
        }
    );

    CLUTCHLOG(progress, "Check data consistency...");
    try {
        tr.check_tables();
        tr.check_genes();
    } catch(const frictionless::DataInconsistent& e) {
        EXIT_ON_ERROR(DataInconsistent, e.what());
    } catch(const frictionless::DataRowFormat& e) {
        EXIT_ON_ERROR(DataRowFormat, e.what());
    }
    CLUTCHLOG(note, "OK -- checks passed.");

}

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

    auto& log = clutchlog::logger();
    ASSERT(log.levels().contains(log_level));
    log.threshold(log_level);
    log.out(std::cerr);
    log.depth_mark(" >");
    log.style(clutchlog::level::critical,
              clutchlog::fmt::fg::black,
              clutchlog::fmt::bg::red,
              clutchlog::fmt::typo::bold);
    log.hfill_style(clutchlog::fmt::fg::black);
    log.hfill_max(100);
    #ifndef NDEBUG
        log.format("{level_letter}:{depth_marks} {msg} {hfill} {func} @ {file}:{line}\n");
    #else
        log.format("{level_letter}:{depth_marks} {msg}\n");
    #endif

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
        frictionless::FriedmanScore frs(tr,2);
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
        CLUTCHLOG(debug, geneset);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Partial eval...");
        frictionless::EvalSwap peval(frs);
        frictionless::Neighbor
            neighbor(geneset.selected.size());
        neighbor.set(1,3);
        peval(geneset, neighbor);
        CLUTCHLOG(debug, neighbor );
    CLUTCHLOG(progress, "Equivalent full eval...");
        neighbor.move(geneset);
        feval(geneset); // Compare with full eval.
        CLUTCHLOG(debug, geneset);
        ASSERT(geneset.fitness() == neighbor.fitness());
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Algo run...");
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

        frictionless::Signature signature(tr.genes_nb());
        signature.select(0);
        signature.select(1);
        feval(signature);
        CLUTCHLOG(progress, "Initial signature: " << signature);

        search(signature);
    CLUTCHLOG(note, "OK");

    CLUTCHLOG(progress, "Done.");
}
