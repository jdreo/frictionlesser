#define CATCH_CONFIG_MAIN

#include <catch2/catch_all.hpp>

#include <frictionless/log.h>

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

    frictionless::clutchlog_config(); // common config

    // Config fos tests:
    auto& log = clutchlog::logger();
    assert(log.levels().contains(log_level));
    log.threshold(log_level);
    // log.hfill_max(100);
    log.strip_calls(5); // Catch2 inserts a lot of intermediate calls.

    return session.run();
}
