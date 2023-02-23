
#include <gtest/gtest.h>

#include <random>

#include "Dijkstra.h"
#include "SpecialGraphs.h"

class DijkstraTest : public ::testing::Test{};

TEST_F(DijkstraTest, EarlyStopping)
{
    constexpr NodeType size = 100;
    Path<size> path;
    Dijkstra<Path<size>> dijkstra(path);

    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_int_distribution<NodeType> node(0, size - 1);

    NodeType u = node(gen);
    NodeType v = node(gen);

    NodeType s = std::min(u,v);
    NodeType t = std::max(u,v);

    dijkstra.compute(s, t);

    ASSERT_EQ(dijkstra.getDistance(), t-s);
}