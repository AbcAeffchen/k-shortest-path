
#include <gtest/gtest.h>

#include "BFSTree.h"
#include "Graph.h"
#include "SpecialGraphs.h"
#include "tools/GraphGenerator.h"
#include "tools/GraphRW.h"

class GraphTest : public ::testing::Test {};

TEST_F(GraphTest, EdgeSizes)
{
    std::cout << "Max. num nodes (dynamic): " << (1ull << (sizeof(NodeType)*8 -1) ) << std::endl;
    std::cout << "Max. num nodes (static): " << (1ull << (sizeof(NodeType)*8) ) << std::endl;

    // Check the sizes of the different edges. Make sure unused members are not in the data structure.
    std::cout << "Directed, Weighted (float): " << sizeof(Edge<float, true>) << std::endl;
    ASSERT_EQ(sizeof(Edge<float, true>), 2 * std::max(sizeof(NodeType), sizeof(float)));
    std::cout << "Directed, Weighted (double): " << sizeof(Edge<double, true>) << std::endl;
    ASSERT_EQ(sizeof(Edge<double, true>), 2 * std::max(sizeof(NodeType), sizeof(double)));

    std::cout << "Undirected, Weighted (float), dynamic: " << sizeof(Edge<float, false>) << std::endl;
    ASSERT_EQ(sizeof(Edge<float, false>), (3-1) * std::max(sizeof(EdgeIdType), sizeof(float)));   // target and weight packs nicely in the is case
    std::cout << "Undirected, Weighted (double), dynamic: " << sizeof(Edge<double, false>) << std::endl;
    ASSERT_EQ(sizeof(Edge<double, false>), 3 * std::max(sizeof(EdgeIdType), sizeof(double)));

    std::cout << "Directed, Unweighted, dynamic: " << sizeof(Edge<void, true>) << std::endl;
    ASSERT_EQ(sizeof(Edge<void, true>), 1 * sizeof(NodeType));
    std::cout << "Undirected, Unweighted, dynamic: " << sizeof(Edge<void, false>) << std::endl;
    ASSERT_EQ(sizeof(Edge<void, false>), 2 * sizeof(EdgeIdType));
}

TEST_F(GraphTest, ConceptValidations)
{
    // Weight concepts
    ASSERT_TRUE(WeightConcept<float>);
    ASSERT_TRUE(WeightConcept<double>);
    ASSERT_FALSE(WeightConcept<int>);
    ASSERT_FALSE(WeightConcept<long>);
    ASSERT_FALSE(WeightConcept<size_t>);
    ASSERT_FALSE(WeightConcept<bool>);
    ASSERT_FALSE(WeightConcept<std::string>);

    // Edge concepts
    ASSERT_FALSE(EdgeConcept<int>);
    using EdgeTypeWD = Edge<float, true>;
    ASSERT_TRUE(EdgeConcept<EdgeTypeWD>);
    using EdgeTypeWU = Edge<float, false>;
    ASSERT_TRUE(EdgeConcept<EdgeTypeWU>);
    using EdgeTypeUD = Edge<void, true>;
    ASSERT_TRUE(EdgeConcept<EdgeTypeUD>);
    using EdgeTypeUU = Edge<void, false>;
    ASSERT_TRUE(EdgeConcept<EdgeTypeUU>);

    // make sure the constructors compile
    Edge<float, true> edge1(0, 1.0);
    Edge<void, true> edge2(0);
    Edge<float, true> edge3(0, 1.0);
    Edge<void, true> edge4(0);
}

/**
 * Generates a path.
 * @param numNodes
 * @return
 */
GraphGenerator::AdjacencyListType generatePath(const NodeType numNodes)
{
    GraphGenerator::AdjacencyListType graph(numNodes);

    for(NodeType u = 0; u < numNodes - 1; u++)
        graph[u].emplace_back(u + 1, 1.0);

    return graph;
}

/**
 * Generates a circle
 * @param size
 * @return
 */
GraphGenerator::AdjacencyListType generateCircle(const NodeType numNodes)
{
    auto graph = generatePath(numNodes);
    graph[numNodes - 1].emplace_back(0, 1.0);

    return graph;
}

void printGraph(const GraphGenerator::AdjacencyListType& graph)
{
    for(size_t u = 0; u < graph.size(); u++)
    {
        std::cout << u << ":";
        for(auto v : graph[u])
            std::cout << " " << v.getTarget();
        std::cout << "\n";
    }
}

TEST_F(GraphTest, RelabelingNodes)
{
    GraphGenerator::AdjacencyListType graph = generatePath(10);

    std::cout << "\nBefore: \n";
    printGraph(graph);

    std::vector<NodeType> removeNodes = {0};
    GraphGenerator::filterGraph(graph, removeNodes);

    std::cout << "\nAfter: \n";
    printGraph(graph);
}

TEST_F(GraphTest, SpecialGraphs)
{
    using BGT = Graph<Edge<void, true>>;
    ASSERT_TRUE(GraphConcept<BGT>);

    using GT = Path<5>;
    ASSERT_TRUE(GraphConcept<GT>);
    GT g;

    GT::NodeType current = 0;
    auto next = g.getOutEdges(current);
    while(!next.empty())
    {
        current = next[0].getTarget();
        next = g.getOutEdges(current);
    }
}

TEST_F(GraphTest, Gilbert)
{
    using GT = Graph<Edge<float, true>, false>;
    unsigned n = 1000;
    unsigned d = 4; // = (n-1) * p (target out-degree)
    double p = static_cast<double>(d) / static_cast<double>(n-1);

    // generates a G(n,p) graph and returns the giant component.
    GT graph = GT::create(GraphGenerator::gilbert(true, true, n, d));

    double numExpectedEdges = graph.getNumNodes() * d;
    double variance = graph.getNumNodes() * p * (graph.getNumNodes() - 1) * (1 - p);
    double sigma = std::sqrt(variance);
    ASSERT_LE(graph.getNumNodes(), n);
    ASSERT_GE(graph.getNumEdges(), numExpectedEdges - 2.1*sigma);
    ASSERT_LE(graph.getNumEdges(), numExpectedEdges + 2.1*sigma);

    size_t maxDeg = 0;
    std::array<size_t, 5> dist = {0,0,0,0,0};

    double numExpectedEdgesPerNode =  d;
    double variance2 = (graph.getNumNodes()-1) * p * (1 - p);
    double sigma2 = std::sqrt(variance2);

    for(NodeType i = 0; i < graph.getNumNodes(); i++)
    {
        auto degi = graph.getOutEdges(i).size();
        maxDeg = std::max(maxDeg, degi);
        auto tmp = std::floor(static_cast<double>(std::abs(static_cast<long>(degi) - static_cast<long>(numExpectedEdgesPerNode))) / sigma2);
        if(tmp < 4)
            dist[static_cast<size_t>(tmp)]++;
        else
            dist[4]++;
    }
    std::cout << "Max Deg: " << maxDeg << std::endl;
    std::cout << "num nodes with at most i sigma (" << sigma2 << ") distance to the expected value (" << d <<"):" << std::endl;
    for(size_t i = 0; i < 5; i++)
    {
        std::cout << "i = " << i + 1 << ": " << dist[i] << " (" << static_cast<double>(dist[i] ) * 100.0 / graph.getNumNodes() << "%)" << std::endl;
    }
    ASSERT_LE(dist[0], static_cast<double>(graph.getNumEdges()) * 0.68);
    ASSERT_LE(dist[0] + dist[1], static_cast<double>(graph.getNumEdges()) * 0.95);
    ASSERT_LE(dist[2] + dist[3] + dist[4], static_cast<double>(graph.getNumEdges()) * 0.5);
}