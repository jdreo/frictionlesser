#include <string>

#include <eo>


int main(int argc, char* argv[])
{
    eoParser parser(argc, argv);
    eoState state;

    const std::string ranks = parser.createParam<std::string>("unknown", "ranks",
        "Filename of the ranks table", 'r', "Data").value();

    const std::string signatures = parser.createParam<std::string>("unknown", "signatures",
        "Filename of the signatures", 'i', "Data").value();

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

    const bool permute = parser.createParam<bool>(false, "permute",
        "Randomly permute the data to get rid of the signal", 'R', "Data").value();


}
