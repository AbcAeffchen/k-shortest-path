/*
 * Check if a Bitfield is faster than a linear search for "removing" nodes from a graph.
 *
 * Result: for 2^25 nodes and less a bitfield is faster. But construction the bitfield takes some time. The time it
 * takes to construct the bitfield is enough to make up for the difference in iterating over the neighbors of a few
 * thousand nodes. for 2^26 nodes and more a bitfield is slower and takes a lot of time to construct (>250.000 ns).
 */

#include <benchmark/benchmark.h>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_set>
#include <utility>

#include "concepts/graph.h"



class DynamicGraphSpanBaseline
{
public:
    NodeType currentDeviationNodeIndex;

public:
    DynamicGraphSpanBaseline(const std::vector<NodeType>& /*forbiddenEdges*/,
                                  const std::vector<NodeType>& /*parentPath*/,
                                  NodeType /*parentDeviationNodeIndex*/,
                                  NodeType currentDeviationNodeIndex,
                                  NodeType /*numNodes*/,
                                  NodeType /*avgDegree*/)
        : //forbiddenEdges(forbiddenEdges),
          //parentPath(parentPath),
          //parentDeviationNodeIndex(parentDeviationNodeIndex),
          currentDeviationNodeIndex(currentDeviationNodeIndex)
          //nodeRemoved(numNodes, false)
    {}

    auto getOutEdges(const NodeType /*u*/, const std::span<const NodeType>& allOutEdges) const
    {
        return allOutEdges;
    }
};

std::random_device rd;
std::default_random_engine gen(rd());

void assignRandomNodes(std::vector<NodeType>& neighbors, std::uniform_int_distribution<NodeType> randomNodes)
{
    for(NodeType& e : neighbors)
        e = randomNodes(gen);
}

template<typename DynamicGraphType>
static void BM_IterateEdges(benchmark::State& state)
{
    // Perform setup
    const auto numNodes = static_cast<NodeType>(state.range(1));
    const auto avgDegree = static_cast<NodeType>(state.range(0));
    const auto pathLength = static_cast<NodeType>(std::log2(numNodes) + 2);

    std::uniform_int_distribution<NodeType> randomNodes(0, numNodes - 1);
    std::uniform_int_distribution<NodeType> randomDeviationNodeIndex(0, std::max(pathLength / 2 - 1, 1u));

    // generate parent path data. Fixed for the whole benchmark
    std::vector<NodeType> forbiddenEdges(avgDegree / 4 + 1);
    assignRandomNodes(forbiddenEdges, randomNodes);
    std::vector<NodeType> path(pathLength);
    assignRandomNodes(path, randomNodes);
    NodeType parentDeviationNodeIndex = randomDeviationNodeIndex(gen);
    NodeType currentDeviationNodeIndex = parentDeviationNodeIndex + randomDeviationNodeIndex(gen);

    // construct dynamic graph object. It is not part of this benchmark.
    DynamicGraphType dynamicGraph(forbiddenEdges, path, parentDeviationNodeIndex, currentDeviationNodeIndex, numNodes, avgDegree);
    std::vector<NodeType> neighbors(avgDegree);
    NodeType tmp;
    // this code is timed
    for (auto _ : state)
    {
        state.PauseTiming();
        assignRandomNodes(neighbors, randomNodes);
        dynamicGraph.currentDeviationNodeIndex = parentDeviationNodeIndex + randomDeviationNodeIndex(gen);
        state.ResumeTiming();
        for(const NodeType& e : dynamicGraph.getOutEdges(randomNodes(gen), neighbors))
        {
            benchmark::DoNotOptimize(tmp = e);
        }
    }
}

template<typename DynamicGraphType>
static void BM_ConstructDynamicGraph(benchmark::State& state)
{
    // Perform setup
    const auto numNodes = static_cast<NodeType>(state.range(1));
    const auto avgDegree = static_cast<NodeType>(state.range(0));
    const auto pathLength = static_cast<NodeType>(std::log2(numNodes) + 2);

    std::uniform_int_distribution<NodeType> randomNodes(0, numNodes - 1);
    std::uniform_int_distribution<NodeType> randomDeviationNodeIndex(0, std::max(pathLength / 2 - 1, 1u));

    // generate parent path data. Fixed for the whole benchmark
    std::vector<NodeType> forbiddenEdges(avgDegree / 4 + 1);
    std::vector<NodeType> path(pathLength);
    NodeType parentDeviationNodeIndex = randomDeviationNodeIndex(gen);
    NodeType currentDeviationNodeIndex = parentDeviationNodeIndex + randomDeviationNodeIndex(gen);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(DynamicGraphType(forbiddenEdges, path, parentDeviationNodeIndex, currentDeviationNodeIndex, numNodes, avgDegree));
    }
}

//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupHashSet)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupBinarySearch)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupLinearSearch<true>)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupLinearSearch<false>)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupLinearSearch1)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupLinearSearch2)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphSpanBaseline)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 26, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphGeneratorBaseline)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 26, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupNodeBitField)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_IterateEdges, DynamicGraphMockupEdgeBitField)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});

//BENCHMARK_TEMPLATE(BM_ConstructDynamicGraph, DynamicGraphMockupLinearSearch)->RangeMultiplier(2)->Ranges({{64, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_ConstructDynamicGraph, DynamicGraphMockupBinarySearch)->RangeMultiplier(2)->Ranges({{64, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_ConstructDynamicGraph, DynamicGraphMockupNodeBitField)->RangeMultiplier(2)->Ranges({{64, 64}, {1 << 20, 1 << 26}});
//BENCHMARK_TEMPLATE(BM_ConstructDynamicGraph, DynamicGraphMockupEdgeBitField)->RangeMultiplier(2)->Ranges({{4, 64}, {1 << 20, 1 << 26}});

// Run the benchmark
BENCHMARK_MAIN();