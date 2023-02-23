
#include <gtest/gtest.h>

#include <random>

#include "tools/GraphRW.h"
#include "DeltaStepping.h"
#include "KSPAlgorithm.h"

class StatsTest : public ::testing::Test{};

template<GraphConcept T>
using SSSPAlgo = DeltaSteppingSingleThread<true, T, SSSPTree>;

template<GraphConcept T>
using SSSPAlgo2w = DeltaSteppingSingleThread<true, T, SSSPTreeTwoWay>;

TEST_F(StatsTest, NumberOfRecords)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/tiny_gnp.metis");

    constexpr unsigned k = 100;

    auto [testGraphGuided, reverseSSSPTree] = testGraph.precomputeSSSPGuiding<SSSPAlgo2w>(0, 1);

    FengsAlgorithm<GraphType, SSSPAlgo, 1, true> feng(testGraphGuided, k);
    feng.computeWithSSSPSkip<false, SSSPTreeTwoWay>(0, reverseSSSPTree.get());

    const auto& stats = feng.getStats();
    const auto& yellowGraphSizes = stats.getYellowGraphSizes();
    const auto& numExpressEdges = stats.getNumExpressEdges();
    const auto& exploredNodes = stats.getExploredNodes();
    const auto totalPathComputations = stats.getTotalPathComputations();

    size_t totalRecordCount = 0;

    ASSERT_EQ(k - 1, yellowGraphSizes.size());
    ASSERT_EQ(k - 1, numExpressEdges.size());
    ASSERT_EQ(k - 1, exploredNodes.size());

    for(unsigned i = 0; i < k - 1; i++)
    {
        totalRecordCount += yellowGraphSizes[i].size();
        ASSERT_EQ(yellowGraphSizes[i].size(), numExpressEdges[i].size()) << i;
        ASSERT_EQ(yellowGraphSizes[i].size(), exploredNodes[i].size()) << i;
    }

    ASSERT_EQ(totalPathComputations, totalRecordCount);

}
