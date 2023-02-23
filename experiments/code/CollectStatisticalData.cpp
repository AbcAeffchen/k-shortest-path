#include <iostream>
#include <random>

#include "CLI11/CLI11.hpp"

#include "BFSTree.h"
#include "tools/GraphRW.h"
#include "tools/Output.h"
#include "tools/KspRunner.h"

int main(int argc, char** argv)
{
    CLI::App app{"CollectStatisticalData"};

    std::string file;
    unsigned numIterations;
    unsigned skipIterations;
    size_t seed;
    double delta;
    unsigned k;
    Algo algo;
    bool useSSSPGuiding;
    unsigned attemptSSSPSkip;
    bool skipByLength;
    bool useEarlyStopping;

    initCli(app, file, numIterations, skipIterations, seed, delta, k,
            algo, useSSSPGuiding, attemptSSSPSkip, skipByLength, useEarlyStopping);

    CLI11_PARSE(app, argc, argv);

    fixCliOptions(algo, useSSSPGuiding, attemptSSSPSkip, skipByLength, seed);

    Print::info() << "Start reading graph from file using delta = " << delta << ". " << std::flush;

    using GraphType = Graph<Edge<double, true>, true>;
    GraphType graph = GraphRW::readMetisFile<double, true, true>(file, delta);

    Print::info() << "Done.\n"
                     "Start calculating k = " << k << " paths\n"
                     "using delta = " << delta
                  << ", skip = " << attemptSSSPSkip
                  << (useSSSPGuiding ? " with" : " without") << " guiding." << std::endl;

    Print::info() << "Seed: " << seed << std::endl;

    std::default_random_engine gen(seed);
    std::uniform_int_distribution<NodeType> drawNode(0, graph.getNumNodes());

    for(unsigned i = 0; i < numIterations; i++)
    {
        NodeType source = drawNode(gen);
        NodeType target;
        do
        {
            target = drawNode(gen);
        }
        while(target == source);

        if(skipIterations > 0)
        {
            skipIterations--;
            Print::info() << "Iteration " << i << " skipped\n";
            continue;
        }
        Print::info() << "Iteration " << i << std::endl;
        try
        {
            KspRunner::run<true, 1, 2>(graph,
                                       algo, useEarlyStopping,
                                       useSSSPGuiding, attemptSSSPSkip, skipByLength,
                                       source, target, k, seed, i);
        }
        catch(...)
        {
            Print::info() << "-> not connected " << std::endl;
        }
        Print::data() << "," << std::endl;
    }

    Print::info() << "\nDone" << std::endl;
}