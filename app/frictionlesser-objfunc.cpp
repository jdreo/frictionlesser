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

    std::string exprsfile;
    std::string ranksfile;
    std::string cachetransfile;
    std::string cachesizefile;
    bool cacheonly;
    double alpha;
    std::string log_level;
    std::string log_file;
    std::string log_func;
    size_t log_depth;
    size_t max_errors;
    double epsilon;

    std::tie(
        exprsfile,
        ranksfile,
        cachetransfile,
        cachesizefile,
        cacheonly,
        alpha,
        log_level,
        log_file,
        log_func,
        log_depth,
        max_errors,
        epsilon
    ) = frictionless::make_parser(argparser);

}
