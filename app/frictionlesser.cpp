#include <string>
#include <iostream>

#include <eo>
#include <mo>

#include <clutchlog/clutchlog.h>

#include <frictionless/frictionless.hpp>


int main(int argc, char* argv[])
{
    eoParser parser(argc, argv);
    // eoState state;

    const std::string ranks = parser.createParam<std::string>("", "ranks",
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

    make_verbose(parser);
    make_help(parser);

    auto& log = clutchlog::logger();
    log.out(std::clog);
    log.depth_mark(">");
    log.threshold(clutchlog::level::xdebug);
    // log.format("[{name}] {level}: {depth_marks} {msg}\n");
    // log.style(clutchlog::level::progress,clutchlog::fmt::fg::blue);

    ASSERT(ranks != "", error);

    std::ifstream ifs;
    ifs.open(ranks);
    ASSERT(ifs.is_open(), error);
    frictionless::RankedTranscriptome rt(ifs);
    ifs.close();

    std::clog << "Loaded a ranked transcriptome of "
              << rt.genes().size() << " genes, "
              << rt.affiliations().size() << " cells, and "
              << rt.samples_nb() << " samples."
              << std::endl;
    //std::cout << rt.ranks_table(true) << std::endl;

    frictionless::Signature null(rt.genes().size(), 0);

    rng.reseed(seed);
    eoUniformGenerator<frictionless::Signature::AtomType> unigen;
    eoInitFixedLength<frictionless::Signature> rinit(rt.genes().size(), unigen);

    frictionless::Signature alea(rt.genes().size(),0);
    rinit(alea);

    alea.printOn(std::cout);
    std::cout << std::endl;

    for(size_t i=0; i < alea.size(); ++i) {
        if(alea[i]) {
            std::cout << rt.gene(i) << " ";
        }
    }
    std::cout << std::endl;
}
