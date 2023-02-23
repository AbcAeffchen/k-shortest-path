
#include <gtest/gtest.h>

#include "DynamicGraph.h"
#include "tools/GraphRW.h"

class DynamicGraphTest : public ::testing::Test {};

TEST_F(DynamicGraphTest, CustomEdgeRange)
{
    auto graph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1.0f);
    KSPPath<float> testPath{4, 0, 0, {}, {0, 1, 3, 7, 9}};

    DynamicGraph dynamicGraph(graph, testPath, 0);

    std::cout << "Edges: ";
    auto test = dynamicGraph.getOutEdges<true>(0);
    for(const auto& e : test)
        std::cout << e.getTarget() << ", " << std::flush;

    std::cout << " -";
}