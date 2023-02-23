#ifndef KSP_KIM_H
#define KSP_KIM_H

#include <climits>

#include "omp.h"
#include "concepts/sssp.h"
#include "Graph.h"
#include "DynamicGraph.h"
#include "KspBasics.h"

template<
    GraphConcept GraphType,
    template<typename> typename SSSPAlgorithmType,
    template<typename> typename DynamicGraphType,
    unsigned numThreads = 1
>
requires SSSPAlgoConcept<SSSPAlgorithmType<GraphType>>
    && std::is_same_v<DynamicGraphType<GraphType>, DynamicGraph<GraphType>>
    && (numThreads > 0)
    && (numThreads >= 1)
class KIM
{
private:
    using DistanceType = typename GraphType::WeightType;
    using ForbiddenEdgesList = typename KSPPath<DistanceType>::ForbiddenEdgesListType;

    const GraphType& graph;
    const unsigned k;
    CandidateList<DistanceType> candidates;
    std::vector<KSPPath<DistanceType>> result;
    std::vector<std::vector<NodeType>> W;

public:
    KIM(const GraphType& graph, const unsigned k) noexcept
      : graph(graph), k(k), candidates(k), W(k)
    {
        result.reserve(k);
    }

public:
    void compute(const NodeType s, const NodeType t) noexcept requires std::is_same_v<DynamicGraphType<GraphType>, DynamicGraph<GraphType>>
    {
        // get shortest path
        SSSPAlgorithmType<GraphType> ssspAlgo(graph);
        ssspAlgo.compute(s, t);

        if(!ssspAlgo.pathFound())
            return;

        result.push_back({ssspAlgo.getDistance(), 0, 0, {}, ssspAlgo.getPath()});
        W[0] = {0, result[0].path.size() - 1};

        if(k < 2)
            return;

        // compute second-shortest path. FSP adds it to the candidate list.
        FSP(graph, s, t, result[0].path, {});

        _compute(t);
    }

private:




    // todo make a parallel version of _compute
    template<template<typename> typename SSSPTreeType = SSSPTree>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    void _compute(const NodeType t) noexcept
    {
        std::unordered_map<uint64_t, std::vector<NodeType> > forbiddenNeighbors;

        for(unsigned i = 1; i < k; i++)
        {
            if(!candidates.hasNext())
                return;

            result << candidates;

            if(i == k - 1)
                return;

            W[i] = {result[i].deviationNodeIndex + 1, result[i].path.size() - 1};

            auto parentIndex = result[i].parentPathId;
            auto deviationNodeIndex = result[i].deviationNodeIndex;
            auto& B = forbiddenNeighbors[hashBIndex(parentIndex, deviationNodeIndex)];
            B.push_back(result[i].path[deviationNodeIndex + 1]);

            // Compute P_a
            if(deviationNodeIndex + 1 < result[i].path.size() - 1)
            {
                DynamicGraphType<GraphType> G_a(graph, result[i], deviationNodeIndex, nullptr);
                FSP(G_a, result[i].path[deviationNodeIndex + 1], t,
                    std::span<NodeType>(result[i].path.begin() + deviationNodeIndex + 1, result[i].path.end()),
                    std::span<NodeType>(result[i].path.begin(), deviationNodeIndex + 1), i);
            }

            // compute P_b
            {
                auto gamma = std::ranges::min_element(W[parentIndex] | std::views::filter([deviationNodeIndex](auto e){ return e > deviationNodeIndex; }));
                assert(deviationNodeIndex < gamma);
                result[parentIndex].forbiddenEdges = B;
                DynamicGraphType<GraphType> G_b(graph, result[parentIndex], deviationNodeIndex, nullptr);
                FSP(G_b, result[parentIndex].path[deviationNodeIndex], t,
                    std::span<NodeType>(result[parentIndex].path.begin() + deviationNodeIndex, result[parentIndex].path.begin() + gamma +  1),
                    std::span<NodeType>(result[parentIndex].path.begin(), deviationNodeIndex), parentIndex);
            }

            // compute P_c
            if(std::ranges::find(W[parentIndex], deviationNodeIndex) != W[parentIndex].end())
            {
                W[parentIndex].push_back(deviationNodeIndex);
                auto& parentDeviationIndex = result[parentIndex].deviationNodeIndex;
                auto gamma = std::ranges::max_element(W[parentIndex] | std::views::filter([parentDeviationIndex, deviationNodeIndex](auto e){ return parentDeviationIndex < e && e < deviationNodeIndex; }));
                assert(gamma < deviationNodeIndex);
                result[parentIndex].forbiddenEdges = forbiddenNeighbors[hashBIndex(parentIndex, gamma)];
                DynamicGraphType<GraphType> G_c(graph, result[parentIndex], gamma, nullptr);
                FSP(G_c, result[parentIndex].path[gamma], t,
                    std::span<NodeType>(result[parentIndex].path.begin() + gamma, result[parentIndex].path.begin() + deviationNodeIndex +  1),
                    std::span<NodeType>(result[parentIndex].path.begin(), gamma), parentIndex);
            }
        }
    }

    void FSP(const DynamicGraphType<GraphType>& tmpGraph, const NodeType s, const NodeType t,
             const std::span<NodeType> prefixPath, const std::span<NodeType> totalPrefixPath,
             const unsigned parentPathIndex) noexcept
    {
        size_t maxDeviationIndex = prefixPath.size();
        SSSPAlgorithmType<DynamicGraphType<GraphType>> ssspForward(tmpGraph);
        ssspForward.comute(s);
        // todo compute tree children and deviation indexes
        SSSPAlgorithmType<DynamicGraphType<GraphType>> ssspBackward(tmpGraph);
        ssspBackward.compute(t);
        // todo compute tree children and deviation indexes

        // the length of the current shortest path
        DistanceType L = Distance<DistanceType>::max;
        // the end of the prefix path and the start of the suffix path of the shortest path.
        std::pair<NodeType, NodeType> H = {nullNode, nullNode};

        // SEP (stack + while-loop)
        std::vector<NodeType> nodeStack = {s};
        while(!nodeStack.empty())
        {
            const auto u = nodeStack.back();
            nodeStack.pop_back();
            const auto& ssspChildren = ssspForward.getChildren(u);

            size_t prefixDeviationIndex = ssspForward.getDeviationIndex(u);
            size_t suffixDeviationIndex = ssspBackward.getDeviationIndex(u);

            for(auto v : ssspChildren)
                if(ssspForward.getDeviationIndex(v) < maxDeviationIndex)
                    nodeStack.push_back(v);

            if(prefixDeviationIndex < suffixDeviationIndex) // Type 1
            {
                const auto l = ssspForward.getDistance(u) + ssspBackward.getDistance(u);
                if(l < L)
                {
                    L = l;
                    H = {u, u};
                }
            }
            else if(prefixDeviationIndex == suffixDeviationIndex) // Type 1
            {
                for(const auto v : tmpGraph.getNeighbors(u))
                {
                    if(ssspChildren.contains(v) || ssspForward.getDeviationIndex(v) >= maxDeviationIndex)
                        continue;

                    const auto l = ssspForward.getDistance(u) + tmpGraph.getEdgeWeight(u, v) + ssspBackward.getDistance(v);
                    if(l < L)
                    {
                        L = l;
                        H = {u, v};
                    }
                }
            }
        }

        if(L == Distance<DistanceType>::max)
            return;

        addCandidatePath(parentPathIndex, ssspForward.getDeviationIndex(H.first),
                         totalPrefixPath + ssspForward.getShortestPath(H.first) + ssspBackward.getReverseShortestPath(H.second));
    }

    void addCandidatePath(const unsigned i, const unsigned deviationNodeIndex, const PathType& candidate) noexcept
    {
        candidates.push({
            graph.getPathLength(candidate),
            i, deviationNodeIndex,
            result[i].deviationNodeIndex == deviationNodeIndex
                ? result[i].getExtendedList(result[i].path[deviationNodeIndex + 1])
                : ForbiddenEdgesList{result[i].path[deviationNodeIndex + 1]},
            candidate
        });
    }

    static uint64_t hashBIndex(unsigned i, NodeType deviationNodeIndex) noexcept
    {
        constexpr unsigned NodeTypeSize = CHAR_BIT * sizeof(NodeType);
        static_assert(NodeTypeSize <= 32);
        return (static_cast<uint64_t>(i) << NodeTypeSize) + static_cast<uint64_t>(deviationNodeIndex);
    }

public:
    [[nodiscard]] const std::vector<KSPPath<DistanceType>>& getResults() const noexcept
    {
        return result;
    }
};

#endif //KSP_KIM_H
