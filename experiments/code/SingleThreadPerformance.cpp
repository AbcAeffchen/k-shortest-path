#include <iostream>

#include "CLI11/CLI11.hpp"

#include "BFSTree.h"
#include "tools/GraphRW.h"
#include "tools/Timer.h"

#include "tools/KspRunner.h"

int main(int argc, char** argv)
{
    CLI::App app{"SingleThreadPerformance"};

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

    initCli<false>(app, file, numIterations, skipIterations, seed, delta, k,
            algo, useSSSPGuiding, attemptSSSPSkip, skipByLength, useEarlyStopping);

    NodeType source, target;
    app.add_option("-s,--source", source, "Source node.")
        ->required()
        ->check(CLI::NonNegativeNumber);
    app.add_option("-t,--target", target, "Target node.")
        ->required()
        ->check(CLI::NonNegativeNumber);

    CLI11_PARSE(app, argc, argv);

    fixCliOptions(algo, useSSSPGuiding, attemptSSSPSkip, skipByLength, seed);

    Print::info() << "Start reading graph from file using delta = " << delta << ". " << std::flush;

    using GraphType = Graph<Edge<double, true>, true>;
    GraphType graph = GraphRW::readMetisFile<double, true, true>(file, delta);
    Timer timer;

    Print::info() << "Done.\n"
                     "Start calculating k = " << k << " paths\n"
                     "using delta = " << delta
                  << ", skip = " << attemptSSSPSkip
                  << (useSSSPGuiding ? " with" : " without") << " guiding." << std::endl;

    {
        ScopedTimer st(timer);

        KspRunner::run<false, 1, 1>(graph,
                                    algo, useEarlyStopping,
                                    useSSSPGuiding, attemptSSSPSkip, skipByLength,
                                    source, target, k);
    }

    Print::info() << "Done. " << timer << std::endl;
}