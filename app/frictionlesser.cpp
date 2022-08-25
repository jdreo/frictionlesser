#include <string>

#include <eo>


int main(int argc, char* argv[])
{
    eoParser parser(argc, argv);
    eoState state;

    std::string ranknames = parser.createParam<std::string>("unknown", "RankNames",
        "Table of rank names.", 'r', "Data").value();
}
