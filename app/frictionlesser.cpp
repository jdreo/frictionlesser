#include <string>
#include <iostream>
#include <fstream>

#include <eo>
#include <mo>

#include <clutchlog/clutchlog.h>

#include <frictionless/frictionless.hpp>


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
            auto& log = clutchlog::logger(); \
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

#define EXIT(err_code) exit(static_cast<unsigned char>(Error::err_code))


int main(int argc, char* argv[])
{
    eoParser parser(argc, argv);
    // eoState state;

    const std::string ranksfile = parser.createParam<std::string>("", "ranks",
        "Filename of the ranks table", 'r', "Data", true).value();

    const std::string signatures = parser.createParam<std::string>("", "signatures",
        "Name of a file containing candidate/starting signatures", 'i', "Data").value();

    const bool permute = parser.createParam<bool>(false, "permute",
        "Randomly permute the data to get rid of the signal", 'R', "Data").value();

    const size_t genes = parser.createParam<size_t>(50, "genes",
        "Number of genes in the signatures", 'g', "Parameters").value();

    const double alpha = parser.createParam<double>(1, "alpha",
        "Score adjustment exponent on the number of genes", 'a', "Parameters").value();

    const double beta = parser.createParam<double>(2, "beta",
        "Exponent on the log of p-values", 'b', "Parameters").value();

    const bool optimum = parser.createParam<bool>(true, "optimum",
        "Stop search only when having reach a local optimum", 'o', "Stopping Criterion").value();

    const long seed = parser.createParam<long>(0, "seed",
        "Seed of the pseudo-random generator (0 = Epoch)", 's', "Misc").value();

    const std::string log_level = parser.createParam<std::string>("Progress", "log-level",
        "Maximum depth level of logging (Critical<Error<Warning<Progress<Note<Info<Debug<XDebug, default=Progress)", 'l', "Misc").value();

    const size_t max_errors = parser.createParam<size_t>(30, "max-errors",
        "Maximum number of errors reported for each check", 'e', "Misc").value();

    make_verbose(parser);
    make_help(parser);

    auto& log = clutchlog::logger();
    ASSERT(log.levels().contains(log_level));
    log.threshold(log_level);
    log.out(std::cerr);
    log.depth_mark(">");
    log.style(clutchlog::level::critical,clutchlog::fmt::fg::black,clutchlog::fmt::bg::red,clutchlog::fmt::typo::bold);

    if(not std::filesystem::exists(ranksfile)) {
        EXIT_ON_ERROR(No_File, "Input ranks data file does not exists."); }

    CLUTCHLOG(progress, "Open `" << ranksfile << "`...");

    std::ifstream ifs(ranksfile);
    if(ifs.fail()) {
        EXIT_ON_ERROR(Unreadable, "Input ranks data file cannot be read."); }
    ASSERT(ifs.is_open());
    frictionless::RankedTranscriptome* rt;
    try {
        rt = new frictionless::RankedTranscriptome(ifs, max_errors);
        ifs.close();
    } catch(const frictionless::DataInconsistent& e) {
        EXIT_ON_ERROR(DataInconsistent, e.what());
    } catch(const frictionless::DataRowFormat& e) {
        EXIT_ON_ERROR(DataRowFormat, e.what());
    } catch(const frictionless::DataSumRanks& e) {
        EXIT_ON_ERROR(DataSumRanks, e.what());
    }
    CLUTCHLOG(progress, "Loaded "
        << rt->genes().size() << " genes, "
        << rt->affiliations().size() << " cells, and "
        << rt->samples_nb() << " samples.");

    CLUTCHCODE(xdebug,
        if( rt->ranks().size() < 90 and rt->genes_nb() < 125) {
            CLUTCHLOG(xdebug, "Loaded ranks table:" << std::endl << rt->format_ranks(true) );
        }
    );
    CLUTCHCODE(debug,
        CLUTCHLOG(debug, "Samples:");
        for(const auto& i : rt->samples()) {
            CLUTCHLOG(debug, "\t- " << rt->sample_name(i) << ": " << rt->cells(i).size() << " cells" );
        }
    );

    CLUTCHLOG(debug, "Test signatures data structures...");
    frictionless::Signature null(rt->genes().size(), 0);

    rng.reseed(seed);
    eoUniformGenerator<frictionless::Signature::AtomType> unigen;
    eoInitFixedLength<frictionless::Signature> rinit(rt->genes().size(), unigen);

    frictionless::Signature alea(rt->genes().size(),0);
    rinit(alea);

    alea.printOn(std::clog);
    std::clog << std::endl;

    for(size_t i=0; i < alea.size(); ++i) {
        if(alea[i]) {
            std::clog << rt->gene_name(i) << " ";
        }
    }
    std::clog << std::endl;
    CLUTCHLOG(debug, "OK");

    CLUTCHLOG(progress, "Pre-compute Friedman score cache...");
    try {
        frictionless::FriedmanScore fs(*rt,2);
    } catch(const frictionless::DataInconsistent& e) {
        EXIT_ON_ERROR(DataInconsistent, e.what());
    } catch(const frictionless::DataRowFormat& e) {
        EXIT_ON_ERROR(DataRowFormat, e.what());
    } catch(const frictionless::DataSumRanks& e) {
        EXIT_ON_ERROR(DataSumRanks, e.what());
    }
    CLUTCHLOG(note, "OK");
}
