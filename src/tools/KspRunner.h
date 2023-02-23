#ifndef SRC_TOOLS_KSPRUNNER
#define SRC_TOOLS_KSPRUNNER

#include <random>

#include "CLI11/CLI11.hpp"

#include "KSPAlgorithm.h"
#include "DeltaStepping.h"

#ifndef THREADS_KSP
#define THREADS_KSP 1
#endif

#ifndef THREADS_DS
#define THREADS_DS 1
#endif

unsigned constexpr numKspThreads = std::max(THREADS_KSP, 1); // > 0 ? THREADS_KSP : 1;
unsigned constexpr numDsThreads = std::max(THREADS_DS, 1); //THREADS_DS > 0 ? THREADS_DS : 1;

static_assert(numKspThreads > 0);
static_assert(numDsThreads > 0);
static_assert(numDsThreads * numKspThreads <= 64);

template<unsigned numThreads, bool useEarlyStopping>
struct DS
{
    template<typename GT>
    using type = DeltaSteppingMultiThread<numThreads, useEarlyStopping, GT, SSSPTree>;
};

template<unsigned numThreads, bool useEarlyStopping>
struct DSTW
{
    template<typename GT>
    using type = DeltaSteppingMultiThread<numThreads, useEarlyStopping, GT, SSSPTreeTwoWay>;
};

enum class Algo{Yen, Feng};

class KspRunner
{
public:
    template<bool collectStatistics, unsigned kspThreads, unsigned dsThreads, GraphConcept GraphType>
    requires (!collectStatistics || kspThreads == 1)
             && (kspThreads > 0) && (dsThreads > 0)
             && (kspThreads * dsThreads <= 64)
    static void run(const GraphType& graph,
             const Algo algo, const bool useEarlyStopping, const bool useSSSPGuiding, const unsigned attemptSSSPSkip, const bool skipByLength,
             const NodeType source, const NodeType target, const unsigned k,
             [[maybe_unused]] const size_t seed = 0, [[maybe_unused]] const unsigned i = 0)
    {
        if(useEarlyStopping)
        {
            if(skipByLength)
                _run<collectStatistics, kspThreads, dsThreads, GraphType, true, true>(graph,
                                                                                      algo, useSSSPGuiding, attemptSSSPSkip,
                                                                                      source, target, k, seed, i);
            else
                _run<collectStatistics, kspThreads, dsThreads, GraphType, true, false>(graph,
                                                                                       algo, useSSSPGuiding, attemptSSSPSkip,
                                                                                       source, target, k, seed, i);
        }
        else
        {
            if(skipByLength)
                _run<collectStatistics, kspThreads, dsThreads, GraphType, false, true>(graph,
                                                                                       algo, useSSSPGuiding, attemptSSSPSkip,
                                                                                       source, target, k, seed, i);
            else
                _run<collectStatistics, kspThreads, dsThreads, GraphType, false, false>(graph,
                                                                                        algo, useSSSPGuiding, attemptSSSPSkip,
                                                                                        source, target, k, seed, i);
        }
    }

private:
    template<bool collectStatistics, unsigned kspThreads, unsigned dsThreads, GraphConcept GraphType,
             bool useEarlyStopping, bool skipByLength>
    requires (!collectStatistics || kspThreads == 1)
             && (kspThreads > 0) && (dsThreads > 0)
             && (kspThreads * dsThreads <= 64)
    static void _run(const GraphType& graph,
             const Algo algo, const bool useSSSPGuiding, const unsigned attemptSSSPSkip,
             const NodeType source, const NodeType target, const unsigned k,
             [[maybe_unused]] const size_t seed, [[maybe_unused]] const unsigned i)
    {
        if(algo == Algo::Yen)
        {
            if(useSSSPGuiding || attemptSSSPSkip > 0)
            {
                auto[guidedGraph, reverseSSSPTree] = graph.template precomputeSSSPGuiding<DS<dsThreads, useEarlyStopping>::template type>(source, target);

                if(useSSSPGuiding)
                {
                    YensAlgorithm<GraphType, DS<dsThreads, useEarlyStopping>::template type, true, kspThreads, collectStatistics> yen(guidedGraph, k);

                    if(attemptSSSPSkip == 1)
                        yen.template computeWithSSSPSkip<skipByLength, SSSPTree>(source, reverseSSSPTree.get());
                    else if(attemptSSSPSkip == 2)
                        yen.template computeWithExtendedSSSPSkip<skipByLength, SSSPTree>(source, reverseSSSPTree.get());
                    else
                        yen.compute(reverseSSSPTree->getReversePath(source));

                    if constexpr(collectStatistics)
                        yen.printStatisticalData(attemptSSSPSkip, skipByLength, seed, i);
                }
                else
                {
                    YensAlgorithm<GraphType, DS<dsThreads, useEarlyStopping>::template type, false, kspThreads, collectStatistics> yen(graph, k);

                    if(attemptSSSPSkip == 2)
                        yen.template computeWithExtendedSSSPSkip<skipByLength, SSSPTree>(source, reverseSSSPTree.get());
                    else
                        yen.template computeWithSSSPSkip<skipByLength, SSSPTree>(source, reverseSSSPTree.get());

                    if constexpr(collectStatistics)
                        yen.printStatisticalData(attemptSSSPSkip, skipByLength, seed, i);
                }
            }
            else
            {
                YensAlgorithm<GraphType, DS<dsThreads, useEarlyStopping>::template type, false, kspThreads, collectStatistics> yen(graph, k);
                yen.compute(source, target);

                if constexpr(collectStatistics)
                    yen.printStatisticalData(attemptSSSPSkip, skipByLength, seed, i);
            }
        }
        else
        {
            auto[guidedGraph, reverseSSSPTree] = graph.template precomputeSSSPGuiding<DSTW<dsThreads, useEarlyStopping>::template type>(source, target);
            FengsAlgorithm<GraphType, DS<dsThreads, useEarlyStopping>::template type, kspThreads, collectStatistics> feng(guidedGraph, k);

            if(attemptSSSPSkip == 2)
                feng.template computeWithExtendedSSSPSkip<skipByLength, SSSPTreeTwoWay>(source, reverseSSSPTree.get());
            else
                feng.template computeWithSSSPSkip<skipByLength, SSSPTreeTwoWay>(source, reverseSSSPTree.get());

            if constexpr(collectStatistics)
                feng.printStatisticalData(attemptSSSPSkip, skipByLength, seed, i);
        }
    }
};

template<bool addOperationalSettings = true>
void initCli(CLI::App& app, std::string& file,
             [[maybe_unused]] unsigned& numIterations, [[maybe_unused]] unsigned& skipIterations, [[maybe_unused]] size_t& seed,
             double& delta, unsigned& k, Algo& algo, bool& useSSSPGuiding, unsigned& attemptSSSPSkip, bool& skipByLength, bool& useEarlyStopping)
{
    app.get_formatter()->column_width(30);
    app.add_option("-g,--graph", file, "Path where the graph is stored. METIS-files only.")
        ->required()
        ->check(CLI::ExistingFile);

    if constexpr(addOperationalSettings)
    {
        app.add_option("-n,--num-iterations", numIterations, "Number of iterations to do.")
            ->required()
            ->check(CLI::PositiveNumber);

        skipIterations = 0;
        app.add_option("--skip-iterations", skipIterations,
                       "Number of iterations to skip up front to resume a computation for a given seed.")
            ->capture_default_str()
            ->check(CLI::NonNegativeNumber);

        seed = 0;
        app.add_option("-s,--seed", seed, "The seed used for the PRNG. 0 (default) will choose a random seed.")
            ->capture_default_str()
            ->check(CLI::NonNegativeNumber);
    }

    delta = 0.01;
    app.add_option("-d,--delta", delta, "Delta parameter used in DeltaStepping.")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);

    k = 10;
    app.add_option("-k", k, "Number of shortest paths to compute.")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);

    algo = Algo::Yen;
    app.add_option("-a,--algo", algo, "Algorithm used. 0=Yen (default), 1=Feng")
        ->capture_default_str();

    useSSSPGuiding = false;
    app.add_flag("--guided,!--no-guided", useSSSPGuiding, "Use SSSP Guiding. Using this flag activates a preprocessing "
                                                          "step to guide the SSSP algorithm in the right direction.")
        ->capture_default_str();

    attemptSSSPSkip = 0;
    app.add_option("--ssspSkip", attemptSSSPSkip, "Attempt to skip SSSP computations if possible shortest or "
                                                  "second shortest deviation pulled from a precomputed SSSP tree.")
        ->check(CLI::IsMember({0, 1, 2}))
        ->capture_default_str();

    skipByLength = true;
    app.add_flag("--skipByLength,!--no-skipByLength", skipByLength, "If skips attempts are made, this option allows to "
                                                                    "use the deviation length to skip SSSP computations.")

        ->capture_default_str();

    useEarlyStopping = true;
    app.add_flag("--earlyStopping,!--no-earlyStopping", useEarlyStopping, "Enables early stopping for SSSP computations."
                                                                          "If turned off, all SSSP computations compute "
                                                                          "the distance to every reachable node.")
        ->capture_default_str();
}

void fixCliOptions(const Algo algo, bool& useSSSPGuiding, unsigned& attemptSSSPSkip, bool& skipByLength, size_t& seed)
{
    if(algo == Algo::Feng)
    {
        useSSSPGuiding = true;

        if(attemptSSSPSkip == 0)
            attemptSSSPSkip = 1;
    }

    if(attemptSSSPSkip == 0)
        skipByLength = false;

    if(seed == 0)
    {
        std::random_device getSeed;
        seed = getSeed();
    }
}

#endif //SRC_TOOLS_KSPRUNNER
