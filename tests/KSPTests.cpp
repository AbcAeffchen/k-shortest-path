
#include <gtest/gtest.h>

#include <random>

#include "tools/GraphRW.h"
#include "Dijkstra.h"
#include "DeltaStepping.h"
#include "KSPAlgorithm.h"

class KSPTest : public ::testing::Test{};

template<GraphConcept T>
using SSSPAlgoWoES = DeltaSteppingSingleThread<false, T, SSSPTree>;

template<GraphConcept T>
using SSSPAlgo = DeltaSteppingSingleThread<true, T, SSSPTree>;

template<GraphConcept T>
using SSSPAlgo2w = DeltaSteppingSingleThread<true, T, SSSPTreeTwoWay>;

template<typename T>
void print(KSPPath<T> path, unsigned i)
{
    std::cout << i << ". " << path.path[0];
    for(auto node : path.path | std::views::drop(1))
        std::cout << " -> " << node;

    std::cout << " (" << path.length << ", " << path.parentPathId << ", " << path.deviationNodeIndex << ")" << std::endl;
}

template<GraphConcept GraphType>
void checkResults(const std::vector<KSPPath<typename GraphType::DistanceType>>& results, const GraphType& graph)
{
    for(unsigned i = 0; i < results.size(); i++)
    {
        print(results[i], i);

        // check for non-decreasing path length. for identical path length path needs to be different.
        if(i > 0)
        {
            ASSERT_FLOAT_EQ(results[i].length, graph.getPathLength(results[i].path)) << "Path length was computed incorrectly";
            ASSERT_LE(results[i - 1].length, results[i].length) << "Paths are in the wrong order";
            if(results[i - 1].length == results[i].length)
            {
                ASSERT_NE(results[i - 1].path, results[i].path) << "Path found multiple times";
            }

            ASSERT_EQ(results[i].path[0], results[0].path[0]) << "Path has different beginning";
            ASSERT_EQ(results[i].path.back(), results[0].path.back()) << "Path has different end";
        }

        // check for loops
        auto path = results[i].path;
        std::ranges::sort(path);

        for(size_t j = 1; j < path.size(); j++)
        {
            ASSERT_NE(path[j-1], path[j]) << "Path is not loopless";
        }
    }
}

TEST_F(KSPTest, FindKPaths)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    YensAlgorithm<GraphType, SSSPAlgo, false> yens(testGraph, 30);

    yens.compute(0, 9);

    checkResults(yens.getResults(), testGraph);
}

TEST_F(KSPTest, UsingPrecomputation)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo>(0, 9);

    YensAlgorithm<GraphType, SSSPAlgo, true> yens(testGraphGuided, 30);

    yens.compute(reverseSSSPTree->getReversePath(0));

    checkResults(yens.getResults(), testGraphGuided);
}

TEST_F(KSPTest, AttemptSSSPSkip_woSDL)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo>(0, 9);

    YensAlgorithm<GraphType, SSSPAlgo, true> yens(testGraphGuided, 30);

    yens.computeWithSSSPSkip<false, SSSPTree>(0, reverseSSSPTree.get());

    checkResults(yens.getResults(), testGraphGuided);
}

TEST_F(KSPTest, AttemptSSSPSkip_wSDL)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo>(0, 9);

    YensAlgorithm<GraphType, SSSPAlgo, true> yens(testGraphGuided, 30);

    yens.computeWithSSSPSkip<true, SSSPTree>(0, reverseSSSPTree.get());

    checkResults(yens.getResults(), testGraphGuided);
}


TEST_F(KSPTest, AttemptSSSPSkip_RandomGrpah)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/tiny_gnp.metis");

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo>(0, 1);

    YensAlgorithm<GraphType, SSSPAlgo, true> yens(testGraphGuided, 100);

    yens.computeWithSSSPSkip<true, SSSPTree>(0, reverseSSSPTree.get());

    checkResults(yens.getResults(), testGraphGuided);
}

TEST_F(KSPTest, Feng_RandomGrpah)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/tiny_gnp.metis");

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo2w>(0, 1);

    FengsAlgorithm<GraphType, SSSPAlgo> fengs(testGraphGuided, 100);

    fengs.computeWithSSSPSkip<true, SSSPTreeTwoWay>(0, reverseSSSPTree.get());

    checkResults(fengs.getResults(), testGraphGuided);
}

TEST_F(KSPTest, Feng_TinyGrpah)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo2w>(0, 9);

    FengsAlgorithm<GraphType, SSSPAlgo> fengs(testGraphGuided, 100);

    fengs.computeWithSSSPSkip<true, SSSPTreeTwoWay>(0, reverseSSSPTree.get());

    checkResults(fengs.getResults(), testGraphGuided);
}

TEST_F(KSPTest, Parallel)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/tiny_gnp.metis");

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo2w>(0, 1);

    FengsAlgorithm<GraphType, SSSPAlgo, 4> fengs(testGraphGuided, 100);

    fengs.computeWithSSSPSkip<true, SSSPTreeTwoWay>(0, reverseSSSPTree.get());

    checkResults(fengs.getResults(), testGraphGuided);
}

TEST_F(KSPTest, CrossCheck)
{
    using GraphType = Graph<Edge<float, true>, true>;
    // file needs to be generated if needed
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/gnp_100k_400k.metis", 0.01f);
    const unsigned k = 50;
    const NodeType s = 0;
    const NodeType t = 1;
    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo2w>(s, t);

    std::cout << "\n--------Yen wo ES---------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgoWoES, false> yens0(testGraph, k);
    yens0.compute(s, t);
    auto yen0Result = yens0.getResults();
    checkResults(yen0Result, testGraph);

    std::cout << "\n--------Yen---------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, false> yens1(testGraph, k);
    yens1.compute(s, t);
    auto yen1Result = yens1.getResults();
    checkResults(yen1Result, testGraph);

    std::cout << "\n--------Yen + guide---------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, true> yens4(testGraphGuided, k);
    yens4.compute(s, t);
    auto yen4Result = yens4.getResults();
    checkResults(yen4Result, testGraphGuided);

    std::cout << "\n--------Yen + skip---------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, false> yens2(testGraph, k);
    yens2.computeWithSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto yen2Result = yens2.getResults();
    checkResults(yen2Result, testGraph);

    std::cout << "\n-------Yen + skip 2--------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, false> yens22(testGraph, k);
    yens22.computeWithExtendedSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto yen22Result = yens22.getResults();
    checkResults(yen22Result, testGraph);

    std::cout << "\n-------Yen + both--------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, true> yens3(testGraphGuided, k);
    yens3.computeWithSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto yen3Result = yens3.getResults();
    checkResults(yen3Result, testGraphGuided);

    std::cout << "\n-------Yen + both 2--------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, true> yens32(testGraphGuided, k);
    yens32.computeWithExtendedSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto yen32Result = yens32.getResults();
    checkResults(yen32Result, testGraphGuided);

    std::cout << "\n-------Yen + both 2 wo length skip--------" << std::endl;

    YensAlgorithm<GraphType, SSSPAlgo, true> yens33(testGraphGuided, k);
    yens33.computeWithExtendedSSSPSkip<false, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto yen33Result = yens33.getResults();
    checkResults(yen33Result, testGraphGuided);

    std::cout << "\n--------Feng--------" << std::endl;

    FengsAlgorithm<GraphType, SSSPAlgo> fengs(testGraphGuided, k);
    fengs.computeWithSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto fengResult = fengs.getResults();
    checkResults(fengResult, testGraphGuided);

    std::cout << "\n--------Feng 2--------" << std::endl;

    FengsAlgorithm<GraphType, SSSPAlgo> fengs2(testGraphGuided, k);
    fengs2.computeWithExtendedSSSPSkip<true, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto feng2Result = fengs2.getResults();
    checkResults(feng2Result, testGraphGuided);

    std::cout << "\n--------Feng 3--------" << std::endl;

    FengsAlgorithm<GraphType, SSSPAlgo> fengs3(testGraphGuided, k);
    fengs3.computeWithExtendedSSSPSkip<false, SSSPTreeTwoWay>(s, reverseSSSPTree.get());
    auto feng3Result = fengs3.getResults();
    checkResults(feng3Result, testGraphGuided);

    std::cout << "\n--------------------" << std::endl;

    for(unsigned i = 0; i < k; i++)
    {
        ASSERT_EQ(yen1Result[i].path, yen0Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen0Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen2Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen2Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen22Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen22Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen3Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen3Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen32Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen32Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen33Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen33Result[i].length;
        ASSERT_EQ(yen1Result[i].path, yen4Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << yen4Result[i].length;
        ASSERT_EQ(yen1Result[i].path, fengResult[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << fengResult[i].length;
        ASSERT_EQ(yen1Result[i].path, feng2Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << feng2Result[i].length;
        ASSERT_EQ(yen1Result[i].path, feng3Result[i].path) << "i = " << i << ", " << yen1Result[i].length << " vs. " << feng3Result[i].length;
    }
}