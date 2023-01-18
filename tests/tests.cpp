#define CATCH_CONFIG_MAIN

#include <catch2/catch_all.hpp>

#include <clutchlog/clutchlog.h>

int main( int argc, char* argv[] ) {
    Catch::Session session; // There must be exactly one instance

    std::string log_level = "Warning";

    // Build a new parser on top of Catch2's
    using namespace Catch::Clara;
    auto cli
        = session.cli()           // Get Catch2's command line parser
        | Opt( log_level, "log-level" ) // bind variable to a new option, with a hint string
          ["-l"]["--log-level"]         // the option names it will respond to
                  ("Log level");        // description string for the help output

    // Now pass the new composite back to Catch2 so it uses that
    session.cli( cli );

    // Let Catch2 (using Clara) parse the command line
    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 ) // Indicates a command line error
        return returnCode;

    auto& log = clutchlog::logger();
    assert(log.levels().contains(log_level));
    log.threshold(log_level);
    log.out(std::cerr);
    log.depth_mark(">");
    log.style(clutchlog::level::critical,
              clutchlog::fmt::fg::black,
              clutchlog::fmt::bg::red,
              clutchlog::fmt::typo::bold);
    log.style(clutchlog::level::error,
              clutchlog::fmt::fg::black,
              clutchlog::fmt::bg::red);
    log.style(clutchlog::level::warning,
              clutchlog::fmt::fg::black,
              clutchlog::fmt::bg::magenta,
              clutchlog::fmt::typo::bold);
    log.style(clutchlog::level::progress,
              clutchlog::fmt::fg::yellow,
              clutchlog::fmt::typo::bold);
    log.style(clutchlog::level::note,
              clutchlog::fmt::fg::green,
              clutchlog::fmt::typo::bold);
    log.style(clutchlog::level::info,
              clutchlog::fmt::fg::green);
    log.style(clutchlog::level::debug,
              clutchlog::fmt::fg::cyan);
    log.style(clutchlog::level::xdebug,
              clutchlog::fmt::fg::blue);
    log.hfill_max(100);
    log.strip_calls(5); // Catch2 inserts a lot of intermediate calls.
    // log.hfill_style(clutchlog::fmt::fg::black);
    #ifndef NDEBUG
        std::ostringstream format;
        clutchlog::fmt reset(clutchlog::fmt::typo::reset);
        clutchlog::fmt discreet(clutchlog::fmt::fg::black);
        clutchlog::fmt bold(clutchlog::fmt::typo::bold);
        format << "{level_fmt}"
               << "{level_letter}:"
               << "{depth_marks} "
               << bold("{msg}")
               << discreet(" {hfill} ")
               << "{level_fmt}{func}"
               << discreet(" @ ")
               << "{level_fmt}{file}"
               << reset << ":"
               << "{level_fmt}{line}"
               << "\n";
        log.format(format.str());
        // log.format("{level_letter}:{depth_marks} {msg} {hfill} {func} @ {file}:{line}\n");
    #else
        log.format("{level_letter}:{depth_marks} {msg}\n");
    #endif

    return session.run();
}
