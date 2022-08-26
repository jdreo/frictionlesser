#include <string>
#include <iostream>

#include <eo>
#include <mo>

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


    assert(ranks != "");

    std::ifstream ifs;
    ifs.open(ranks);
    assert(ifs.is_open());
    frictionless::RankedTranscriptome rt(ifs);
    ifs.close();

    //std::cout << rt.ranks_table(true) << std::endl;

    frictionless::Signature null(rt.genes().size(), 0);

    rng.reseed(seed);
    eoUniformGenerator<char> unigen;
    eoInitFixedLength<frictionless::Signature> rinit(rt.genes().size(), unigen);

    frictionless::Signature alea;
    rinit(alea);
}
