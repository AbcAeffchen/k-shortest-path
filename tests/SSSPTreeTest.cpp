
#include <gtest/gtest.h>

#include <random>

#include "DeltaStepping.h"
#include "SsspTree.h"
#include "tools/GraphRW.h"

class SSSPTreeTest : public ::testing::Test{};

TEST_F(SSSPTreeTest, NumericTest)
{
    using FloatType = double;

    for(size_t t = 2; t <= 64; t *= 2)
    {
        for(double x = 0.8; x < 1.0; x += 0.001)
        {
            for(size_t i = 20; i <= 30; i++)
            {
                auto n = static_cast<size_t>((1 << i) * x);
                auto nt = (n - 1) / t;
                FloatType clusterFrac = static_cast<FloatType>(t) / static_cast<FloatType>(n);
                auto v1 = static_cast<FloatType>(nt) * clusterFrac;
                ASSERT_EQ(0, static_cast<size_t>(v1))
                                    << i << " -> " << n << ", " << t << ", " << x << ", " << v1;
                auto v2 = static_cast<FloatType>(n - 1) * clusterFrac;
                ASSERT_EQ(t - 1, static_cast<size_t>(v2))
                                    << i << " -> " << n << ", " << t << ", " << x << ", " << v2;
            }
        }
    }
}

TEST_F(SSSPTreeTest, ChildrenConnection)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 0.01f);

    DeltaStepping<GraphType, 2, true, SSSPTreeTwoWay> deltaStepping(testGraph);
    deltaStepping.compute(0);

    auto ssspTree = deltaStepping.getSSSPTree();

    for(NodeType node = 0; node < testGraph.getNumNodes(); node++)
    {
        for(auto child : ssspTree->getChildren(node))
        {
            ASSERT_EQ(node, ssspTree->getPredecessor(child)) << node << ", " << child;
        }
    }
}

TEST_F(SSSPTreeTest, ParallelChildren)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 0.01f);

    DeltaStepping<GraphType, 1, true, SSSPTreeTwoWay> deltaStepping1(testGraph);
    deltaStepping1.compute(0);

    DeltaStepping<GraphType, 4, true, SSSPTreeTwoWay> deltaStepping4(testGraph);
    deltaStepping4.compute(0);

    auto ssspTree1 = deltaStepping1.getSSSPTree();
    auto ssspTree4 = deltaStepping4.getSSSPTree();

    for(NodeType node = 0; node < testGraph.getNumNodes(); node++)
    {
        auto children1 = ssspTree1->getChildren(node);
        auto children4 = ssspTree4->getChildren(node);

        std::ranges::sort(children1);
        std::ranges::sort(children4);

        ASSERT_EQ(children1, children4);
    }
}
