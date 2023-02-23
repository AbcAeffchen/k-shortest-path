
#include <gtest/gtest.h>

#include <random>

#include "DeltaStepping.h"
#include "tools/GraphRW.h"

class DeltaSteppingTest : public ::testing::Test{};

template<unsigned numThreads>
void testDeltaSteppingFullTree()
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    DeltaStepping<GraphType, numThreads, true> deltaStepping(testGraph);
    deltaStepping.compute(0);

    ASSERT_TRUE(deltaStepping.pathFound(9));
    ASSERT_FLOAT_EQ(0, deltaStepping.getDistance(0));
    ASSERT_FLOAT_EQ(1, deltaStepping.getDistance(1));
    ASSERT_FLOAT_EQ(2, deltaStepping.getDistance(2));
    ASSERT_FLOAT_EQ(2, deltaStepping.getDistance(3));
    ASSERT_FLOAT_EQ(3, deltaStepping.getDistance(4));
    ASSERT_FLOAT_EQ(4, deltaStepping.getDistance(5));
    ASSERT_FLOAT_EQ(3, deltaStepping.getDistance(6));
    ASSERT_FLOAT_EQ(3, deltaStepping.getDistance(7));
    ASSERT_FLOAT_EQ(4, deltaStepping.getDistance(8));
    ASSERT_FLOAT_EQ(4, deltaStepping.getDistance(9));
}

TEST_F(DeltaSteppingTest, SingleThread)
{
    testDeltaSteppingFullTree<1>();
}

TEST_F(DeltaSteppingTest, MultiThread2)
{
    testDeltaSteppingFullTree<2>();
}

TEST_F(DeltaSteppingTest, MultiThread4)
{
    testDeltaSteppingFullTree<4>();
}

TEST_F(DeltaSteppingTest, MultiThread8)
{
    testDeltaSteppingFullTree<8>();
}

TEST_F(DeltaSteppingTest, EarlyStopping)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    DeltaStepping<GraphType, 1, true> deltaStepping(testGraph);
    deltaStepping.compute(0, 9);

    ASSERT_TRUE(deltaStepping.pathFound(9));
    ASSERT_FLOAT_EQ(4, deltaStepping.getDistance(9));
}

TEST_F(DeltaSteppingTest, EarlyStopping_disabled)
{
    using GraphType = Graph<Edge<float, true>, true>;
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/directed-graph-small.metis", 1);

    DeltaStepping<GraphType, 1, false> deltaStepping(testGraph);
    deltaStepping.compute(0, 9);

    ASSERT_TRUE(deltaStepping.pathFound(9));
    ASSERT_FLOAT_EQ(4, deltaStepping.getDistance(9));
}

TEST_F(DeltaSteppingTest, CheckForEqualResults)
{
    using GraphType = Graph<Edge<float, true>, true>;
    // file needs to be generated if needed
    GraphType testGraph = GraphRW::readMetisFile<float, true, true>("../../tests/data/gnp_1m_4m.metis", 0.01f);

    DeltaStepping<GraphType, 1, true> deltaStepping1(testGraph);
    deltaStepping1.compute(0);

    DeltaStepping<GraphType, 4, true> deltaStepping2(testGraph);
    deltaStepping2.compute(0);

    for(NodeType node = 0; node < testGraph.getNumNodes(); node++)
    {
        ASSERT_FLOAT_EQ(deltaStepping1.getDistance(node), deltaStepping2.getDistance(node));
    }
}