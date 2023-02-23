#include "CLI11/CLI11.hpp"

#include "tools/GraphGenerator.h"
#include "tools/GraphRW.h"
#include "tools/Output.h"

int main (int argc, char *argv[])
{
    CLI::App app{"Graph Generator"};
    app.get_formatter()->column_width(30);
    std::string file;
    app.add_option("-o,--output", file, "Path where the graph is stored")
        ->required()
        ->check(!CLI::ExistingFile);
    bool unweighted = false;
    app.add_flag("--unweighted", unweighted, "If set, the generated graph is unweighted (weighted by default)\n"
                                             "If it is weighted, weights are drawn uniformly at random from [0,1)");
//    bool undirected = false;
//    auto undirectedOption = app.add_flag("--undirected", undirected, "If set, the generated graph is undirected (directed by default)");
    unsigned long seed = 0;
    app.add_option("-s,--seed", seed, "The seed for the random number generator.\n"
                                      "If not provided or zero, a random seed is used.")
        ->check(CLI::NonNegativeNumber);

    NodeType numNodes;

    // G(n,p) Options
    auto gnp = app.add_subcommand("gnp", "Generates a graph following the G(n,p)-Model.");
    gnp->add_option("-n", numNodes, "Number of nodes")->required();
    double avgDegree = 4.0;
    gnp->add_option("-a", avgDegree, "Expected average degree. This translates to p = average degree / number of nodes.\n"
                                     "If the graph is directed, this is the expected out-degree.")
                                     ->capture_default_str()
                                     ->required();
    double directness = 1;
    gnp->add_option("-d,--directness", directness, "The probability for an edge (u,v) that also edge (v,u) is generated.\n"
                                                   "For d = 0 this is a classic undirected graph, for d = 1 it is a classic directed graph.")
           ->capture_default_str()
           ->check(CLI::Range(0.0, 1.0));

    // G(n,m) Options
    auto gnm = app.add_subcommand("gnm", "Generates a graph following the G(n,m)-Model.");
    gnm->add_option("-n", numNodes, "Number of nodes")->required();
    EdgeIdType numEdges;
    gnm->add_option("-m", numEdges, "Number of edges")->required();

    // Grid options
    auto grid = app.add_subcommand("grid", "Generates a 3D grid graph.");
    std::vector<NodeType> dims;
    grid->add_option("-d", dims, "X, Y dimensions of the grid in number of nodes.\n"
                                 "One can set a dimension to one to get a one or two dimensional grid\n"
                                 "where a one dimensional grid would be a path and a two dimensional\n"
                                 "grid would be something like a road network.")
                                 ->expected(2)
                                 ->required();
    double p = 0.95;
    grid->add_option("-p", p, "Probability for an edge to exist.")
        ->capture_default_str();
    bool gridDirected = true;
    grid->add_flag("--directed,!--undirected", gridDirected, "Sets if the generated grid is directed or undirected.")
        ->capture_default_str();
//    bool wrapX, wrapY, wrapZ;
//    grid->add_flag("-x,--warp-around-x", wrapX, "If set, the first and last nodes in X direction are connected.");
//    grid->add_flag("-y,--warp-around-y", wrapY, "If set, the first and last nodes in Y direction are connected.");
//    grid->add_flag("-z,--warp-around-z", wrapZ, "If set, the first and last nodes in Z direction are connected.");

    app.require_subcommand(1);

    CLI11_PARSE(app, argc, argv)

    if(seed == 0)
    {
        std::random_device getSeed;
        seed = getSeed();
    }

    bool directed = true;

    GraphGenerator::AdjacencyListType graph;
    if(app.got_subcommand(gnp))
    {
        directed = directness > 0.0;

        Print::info() << "Start generating G(n,p) graph with n = " << numNodes << ", p = " << avgDegree << " / " << numNodes << ", and directness of " << directness << std::endl;
        graph = GraphGenerator::gilbert(!unweighted, directness, numNodes, avgDegree, seed);

        numEdges = 0;
        for(const auto& neighbors : graph)
            numEdges += neighbors.size();

        Print::info() << "Finished generating " << numEdges << " edges." << std::endl;
        Print::info() << "Giant component has " << graph.size() << " (" << (static_cast<float>(graph.size()) * 100.0f / static_cast<float>(numNodes) ) << "%) nodes." << std::endl;
    }
    else if(app.got_subcommand(gnm))
    {
        Print::info() << "G(n,m) is not ready yet." << std::endl;   // todo
        return 0;
    }
    else if(app.got_subcommand(grid))
    {
        directed = gridDirected;
        graph = GraphGenerator::grid(!unweighted, gridDirected, dims[0], dims[1], p, seed);

        numEdges = 0;
        for(const auto& neighbors : graph)
            numEdges += neighbors.size();

        Print::info() << "Finished generating " << dims[0] << " x " << dims[1] << " grid with " << numEdges << " edges." << std::endl;
        Print::info() << "Giant component has " << graph.size() << " (" << (static_cast<float>(graph.size()) * 100.0f / static_cast<float>(dims[0] * dims[1])) << "%) nodes." << std::endl;
    }

    Print::info() << "Starting to write graph to " << file << std::endl;
    GraphRW::writeMetisFile(directed, !unweighted, static_cast<NodeType>(graph.size()), numEdges, graph, file, seed);
    Print::info() << "Finished process." << std::endl;

    return 0;
}