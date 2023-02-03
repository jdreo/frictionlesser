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

