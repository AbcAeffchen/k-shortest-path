
#include <gtest/gtest.h>
#include <filesystem>

#include "tools/GraphRW.h"

class GraphRWTest : public ::testing::Test {};

TEST_F(GraphRWTest, LoadDirectedGraphFromMetisFile)
{
    auto graph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 2);

    ASSERT_EQ(graph.getNumNodes(), 10);
    ASSERT_TRUE(graph.isDirected);
    ASSERT_TRUE(graph.isWeighted);

    const std::vector<std::vector<NodeType>> neighbors = {{1, 9, 2},
                                                          {3, 7, 4, 2},
                                                          {1, 5, 8, 6},
                                                          {7},
                                                          {7},
                                                          {8},
                                                          {8},
                                                          {9, 1, 8},
                                                          {9, 7, 2}};
    const std::vector<std::vector<NodeType>> lightNeighbors = {{1},
                                                               {3},
                                                               {6},
                                                               {7},
                                                               {},
                                                               {},
                                                               {8},
                                                               {9},
                                                               {9}};

    for(NodeType i = 0; i < 10; i++)
    {
        for(const auto& e : graph.getOutEdges(i))
        {
            ASSERT_TRUE(neighbors[i].end() != std::ranges::find(neighbors[i], e.getTarget()))
                << "Edge not found: (" << i << ", " << e.getTarget() << ")";
        }
    }

    for(NodeType i = 0; i < 10; i++)
    {
        for(const auto& e : graph.getLightOutEdges(i))
        {
            ASSERT_TRUE(lightNeighbors[i].end() != std::ranges::find(lightNeighbors[i], e.getTarget()))
                << "Edge is not light: (" << i << ", " << e.getTarget() << ")";
        }
    }

    for(NodeType i = 0; i < 10; i++)
    {
        for(const auto& e : graph.getHeavyOutEdges(i))
        {
            ASSERT_TRUE(lightNeighbors[i].end() == std::ranges::find(lightNeighbors[i], e.getTarget()))
                << "Edge is not heavy: (" << i << ", " << e.getTarget() << ")";
        }
    }

    PathType path1 = {0, 1, 3, 7, 9};
    ASSERT_FLOAT_EQ(4.0, graph.getPathLength(path1));
    PathType path2 = {0, 2, 6, 8, 9};
    ASSERT_FLOAT_EQ(5.0, graph.getPathLength(path2));
    PathType path3 = {0, 9};
    ASSERT_FLOAT_EQ(8.0, graph.getPathLength(path3));


}
