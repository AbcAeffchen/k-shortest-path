#include "CLI11/CLI11.hpp"

#include "tools/GraphRW.h"
#include "tools/Output.h"


int main (int argc, char *argv[])
{
    CLI::App app{"Graph Converter"};
    app.get_formatter()->column_width(30);
    std::string inputFile;
    app.add_option("-i,--input", inputFile, "Path to the edgelist file")
        ->required()
        ->check(CLI::ExistingFile);
    std::string outputFile;
    app.add_option("-o,--output", outputFile, "Path where the metis file will be stored")
        ->required()
        ->check(!CLI::ExistingFile);
    std::map<std::string, GraphRW::Weights> map{{"Uniform",    GraphRW::Weights::Uniform},
                                                {"Random",     GraphRW::Weights::Random},
                                                {"SumDeg",     GraphRW::Weights::SumDeg},
                                                {"ProdDeg",    GraphRW::Weights::ProdDeg},
                                                {"InvSumDeg",  GraphRW::Weights::InvSumDeg},
                                                {"InvProdDeg", GraphRW::Weights::InvProdDeg}};
    GraphRW::Weights weights(GraphRW::Weights::Random);
    app.add_option("-w,--weights", weights, "Weight assignment")
        ->required()
        ->transform(CLI::CheckedTransformer(map, CLI::ignore_case));

    bool makeUndirected = false;
    app.add_flag("-u,--make-undirected", makeUndirected, "Read graph as undirected.")
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv)

    Print::info() << "Source file: " << inputFile << std::endl;
    Print::info() << "Target file: " << outputFile << std::endl;
    Print::info() << "Weight type: " << static_cast<int>(weights) << std::endl;
    if(makeUndirected)
        Print::info() << "Read graph as undirected." << std::endl;
    else
        Print::info() << "Read graph as is." << std::endl;

    Print::info() << "Reading file and generating weights..." << std::endl;
    auto adjacencyList = GraphRW::convertEdgelistFile(inputFile, weights);

    Print::info() << "Counting edges..." << std::endl;
    size_t numEdges = 0;
    for(const auto& outEdges : adjacencyList)
    {
        numEdges += outEdges.size();
    }

    Print::info() << "Starting to write graph to " << outputFile << std::endl;
    GraphRW::writeMetisFile(true, true, static_cast<NodeType>(adjacencyList.size()), numEdges, adjacencyList, outputFile);
    Print::info() << "Finished process." << std::endl;

    return 0;
}