#include <iostream>

#include "CLI11/CLI11.hpp"

#include "BFSTree.h"
#include "tools/GraphRW.h"
#include "tools/Output.h"
#include "tools/Timer.h"

#include "tools/KspRunner.h"

int main(int argc, char** argv)
{
    CLI::App app{"Performance Measurements"};

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

//    if(numDsThreads * numKspThreads > 64)
//    {
//        Print::info() << "Total number of requested threads is too big (" << (numDsThreads * numKspThreads) << ")" << std::endl;
//        exit(1);
//    }

    Print::info() << "Start reading graph from file using delta = " << delta << ". " << std::flush;

    using GraphType = Graph<Edge<double, true>, true>;
    const GraphType graph = GraphRW::readFile<double, true, true>(file, delta);

    Print::info() << "Done.\n"
                     "Start calculating k = " << k << " paths\n"
                     "using delta = " << delta
                  << ", skip = " << attemptSSSPSkip
                  << (useSSSPGuiding ? " with" : " without") << " guiding." << std::endl;

    if(seed == 0)
    {
        std::random_device getSeed;
        seed = getSeed();
    }

    Print::info() << "Seed: " << seed << std::endl;

    std::default_random_engine gen(seed);
    std::uniform_int_distribution<NodeType> drawNode(0, graph.getNumNodes() - 1);

    auto filePathEnd = file.find_last_of('/');
    auto fileExtensionEnd = file.find_last_of('.');
    auto fileNameStart = filePathEnd == std::string::npos ? 0 : filePathEnd + 1;
    auto fileNameLength = (fileExtensionEnd == std::string::npos ? file.size() : fileExtensionEnd) - fileNameStart;
    auto fileName = file.substr(filePathEnd == std::string::npos ? 0 : fileNameStart, fileNameLength);

    for(unsigned run = 0; run < numIterations; run++)
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
            Print::info() << "Iteration " << run << " skipped\n";
            continue;
        }

        Print::info() << "Iteration " << run << std::endl;

        JSON result;

        result.add("file", fileName);
        result.add("algo", static_cast<std::string>(algo == Algo::Yen ? "yen" : "feng"));
        result.add("earlyStopping", useEarlyStopping);
        result.add("guided", useSSSPGuiding);
        result.add("skipping", attemptSSSPSkip);
        result.add("skipByLength", skipByLength);
        result.add("numNodes", graph.getNumNodes());
        result.add("numEdges", graph.getNumEdges());
        result.add("seed", seed);
        result.add("run", run);
        result.add("source", source);
        result.add("target", target);
        result.add("k", k);
        result.add("delta", delta);
        result.add("kspThreads", numKspThreads);
        result.add("dsThreads", numDsThreads);

        Timer timer;

        {
            ScopedTimer sc(timer);

            KspRunner::run<false, numKspThreads, numDsThreads>(graph,
                                                               algo, useEarlyStopping,
                                                               useSSSPGuiding, attemptSSSPSkip, skipByLength,
                                                               source, target, k);
        }

        result.add("time", timer.getTotalDurationS());

        Print::data() << result << "," << std::endl;
    }
}