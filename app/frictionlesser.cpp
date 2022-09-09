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
        "Filename of the expressions table to be ranked", 'r', "Data").value();

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

    const long seed = argparser.createParam<long>(0, "seed",
        "Seed of the pseudo-random generator (0 = Epoch)", 's', "Misc").value();

    const std::string log_level = argparser.createParam<std::string>("Progress", "log-level",
        "Maximum depth level of logging (Critical<Error<Warning<Progress<Note<Info<Debug<XDebug, default=Progress)", 'l', "Misc").value();

    const size_t max_errors = argparser.createParam<size_t>(30, "max-errors",
        "Maximum number of errors reported for each check", 'e', "Misc").value();

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
        frictionless::Transcriptome ranked = frictionless::rank(tr, /* print_progress */true);
        CLUTCHLOG(note, "OK");

        CLUTCHCODE(xdebug,
            if( ranked.ranks().size() < 90 and ranked.genes_nb() < 125) {
                CLUTCHLOG(xdebug, "Ranked table:" << std::endl << ranked.as_art(true) );
            }
        );

        CLUTCHLOG(progress, "Check ranks consistency...");
        try {
            ranked.check_ranks();
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
        tr.check_ranks();
    } catch(const frictionless::DataSumRanks& e) {
        EXIT_ON_ERROR(DataSumRanks, e.what());
    }
    CLUTCHLOG(note, "OK -- checks passed.");

    // /*
    CLUTCHLOG(debug, "Test signatures data structures...");
    frictionless::Signature null(tr.genes().size(), 0);

    rng.reseed(seed);
    eoUniformGenerator<frictionless::Signature::AtomType> unigen;
    eoInitFixedLength<frictionless::Signature> rinit(tr.genes().size(), unigen);

    frictionless::Signature alea(tr.genes().size(),0);
    rinit(alea);

    CLUTCHCODE(xdebug,
        alea.printOn(std::clog);
        std::clog << std::endl;

        for(size_t i=0; i < alea.size(); ++i) {
            if(alea[i]) {
                std::clog << tr.gene_name(i) << " ";
            }
        }
        std::clog << std::endl;
    );

    CLUTCHLOG(debug, "OK");
    // */
    CLUTCHLOG(progress, "Pre-compute Friedman score cache...");
    // try {
        frictionless::FriedmanScore fs(tr,2);
        std::cout << fs.score(alea) << std::endl;
    // } catch(...) {
        // EXIT_ON_ERROR(DataInconsistent, e.what());
    // }
    CLUTCHLOG(note, "OK");
    CLUTCHLOG(progress, "Done.");
}
