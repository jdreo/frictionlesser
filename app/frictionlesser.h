#pragma once

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
            CLUTCHLOG(critical, msg); \
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

namespace frictionless {

static void check_consistency(const frictionless::Transcriptome& tr)
{
    CLUTCHCODE(debug,
        CLUTCHLOG(debug, "Samples:");
        for(const auto& i : tr.samples()) {
            CLUTCHLOGD(debug, tr.sample_name(i) << ": " << tr.cells(i).size() << " cells" , 1);
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

std::tuple<
    std::string,
    std::string,
    std::string,
    std::string,
    bool,
    double,
    std::string,
    std::string,
    std::string,
    size_t,
    size_t,
    double
> make_parser(eoParser& argparser)
{

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

    const double alpha = argparser.createParam<double>(1, "alpha",
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
}

} // frictionless
