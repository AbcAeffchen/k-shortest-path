#ifndef KSP_YENSALGORITHM_H
#define KSP_YENSALGORITHM_H

#include "omp.h"
#include "concepts/sssp.h"
#include "Graph.h"
#include "DynamicGraph.h"
#include "KspBasics.h"
#include "YellowGraph.h"
#include "tools/Statistics.h"

template<
    GraphConcept GraphType,
    template<typename> typename SSSPAlgorithmType,
    template<typename> typename DynamicGraphType,
    bool guided,
    unsigned numThreads = 1,
    bool collectStatistics = false
>
requires SSSPAlgoConcept<SSSPAlgorithmType<GraphType>>
    && (std::is_same_v<DynamicGraphType<GraphType>, YellowGraph<GraphType>>
        || std::is_same_v<DynamicGraphType<GraphType>, DynamicGraph<GraphType>>)
    && (numThreads > 0)
    && (!collectStatistics || numThreads == 1)
class KSPAlgorithm
{
private:
    using DistanceType = typename GraphType::WeightType;
    using ForbiddenEdgesList = typename KSPPath<DistanceType>::ForbiddenEdgesListType;

    static constexpr bool usesYellowGraph = std::is_same_v<DynamicGraphType<GraphType>, YellowGraph<GraphType>>;

    const GraphType& graph;
    const unsigned k;
    CandidateList<DistanceType> candidates;
    std::vector<KSPPath<DistanceType>> result;
    Statistics<collectStatistics> stats;

public:
    KSPAlgorithm(const GraphType& graph, const unsigned k) noexcept
      : graph(graph), k(k), candidates(k), stats(k)
    {
        result.reserve(k);
    }

    /**
     * Tries to skip sssp calls by using candidate paths from the reverse shortest path tree.
     */
    template<bool skipByLength, template<typename> typename SSSPTreeType>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    void computeWithSSSPSkip(const NodeType s, SSSPTreeType<DistanceType> const * const reverseSSSPTree) noexcept
    {
        _computeWithSSSPSkip<1, skipByLength, SSSPTreeType>(s, reverseSSSPTree);
    }

    /**
     * Tries to skip sssp calls by using candidate paths from the reverse shortest path tree,
     * also considering the second shortest path, if the first one is not simple
     */
    template<bool skipByLength, template<typename> typename SSSPTreeType>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    void computeWithExtendedSSSPSkip(const NodeType s, SSSPTreeType<DistanceType> const * const reverseSSSPTree) noexcept
    {
        _computeWithSSSPSkip<2, skipByLength, SSSPTreeType>(s, reverseSSSPTree);
    }

private:
    template<unsigned attemptSSSPSkip, bool skipByLength, template<typename> typename SSSPTreeType>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    void _computeWithSSSPSkip(const NodeType s, SSSPTreeType<DistanceType> const * const reverseSSSPTree) noexcept
    {
        PathType shortestPath = reverseSSSPTree->getReversePath(s);

        candidates.push({graph.getPathLength(shortestPath), 0, 0, {}, shortestPath});
        _compute<attemptSSSPSkip, skipByLength, SSSPTreeType>(shortestPath.back(), reverseSSSPTree);
    }

public:

    void compute(const PathType& shortestPath) noexcept requires std::is_same_v<DynamicGraphType<GraphType>, DynamicGraph<GraphType>>
    {
        candidates.push({graph.getPathLength(shortestPath), 0, 0, {}, shortestPath});
        _compute<0, false>(shortestPath.back());
    }

    void compute(const NodeType s, const NodeType t) noexcept requires std::is_same_v<DynamicGraphType<GraphType>, DynamicGraph<GraphType>>
    {
        SSSPAlgorithmType<GraphType> ssspAlgo(graph);
        ssspAlgo.compute(s, t);

        if(!ssspAlgo.pathFound())
            return;

        candidates.push({ssspAlgo.getDistance(), 0, 0, {}, ssspAlgo.getPath()});

        _compute<0, false>(t);
    }

private:
    template<unsigned attemptSSSPSkip, bool skipByLength, template<typename> typename SSSPTreeType = SSSPTree>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    void _compute(const NodeType t, SSSPTreeType<DistanceType> const * const reverseSSSPTree = nullptr) noexcept
    {
        // make sure the algorithm is not called multiple times.
        assert(result.empty());

        if constexpr(attemptSSSPSkip > 0)
        {
            assert(reverseSSSPTree != nullptr);
        }

        result << candidates;

        for(unsigned i = 0; i < k - 1; i++)
        {
            computeDeviations<attemptSSSPSkip, skipByLength, SSSPTreeType>(t, reverseSSSPTree, i);

            if(!candidates.hasNext()) [[unlikely]]
                break;

            result << candidates;
        }
    }

    template<unsigned attemptSSSPSkip, bool skipByLength, template<typename> typename SSSPTreeType = SSSPTree>
    void computeDeviations(const NodeType t,
                           SSSPTreeType<DistanceType> const * const reverseSSSPTree,
                           const unsigned i) noexcept
    requires (numThreads == 1)
    {
        DynamicGraphType<GraphType> tmpGraph(graph, result[i], reverseSSSPTree);

        stats.countTotalPathComputations(result[i].path.size() - result[i].deviationNodeIndex - 1);

        for(unsigned deviationNodeIndex = result[i].deviationNodeIndex; deviationNodeIndex < result[i].path.size() - 1; deviationNodeIndex++)
        {
            tmpGraph.advanceDeviationNodeIndex(deviationNodeIndex);

            if constexpr(usesYellowGraph && collectStatistics)
                stats.storeYellowGraphSize(i, tmpGraph.getNumYellowNodes());

            if constexpr(attemptSSSPSkip > 0)
            {
                if(makeSSSPSkipAttempt<SSSPTreeType, attemptSSSPSkip >= 2, skipByLength>(i, deviationNodeIndex, tmpGraph, *reverseSSSPTree))
                {
                    stats.storeNumExploredNodes(i, 0);
                    if constexpr(usesYellowGraph && collectStatistics)
                        stats.storeNumExpressEdges(i, 0);

                    continue;
                }
            }
            else
            {
                // if skips are not attempted by guessing the shortest deviation at least check if the deviation node
                // has any out edges left and skip SSSP construction and initialisation if so.
                auto outEdges = tmpGraph.template getOutEdges<true>(result[i].path[deviationNodeIndex]);

                // no valid out edge exists => go to the next deviation node
                if(outEdges.begin() == outEdges.end())
                {
                    stats.storeNumExploredNodes(i, 1);
                    continue;
                }
            }

            computeDeviationSuffix(tmpGraph, i, deviationNodeIndex, t);
        }
    }

    template<unsigned attemptSSSPSkip, bool skipByLength, template<typename> typename SSSPTreeType = SSSPTree>
    void computeDeviations(const NodeType t,
                           SSSPTreeType<DistanceType> const * const reverseSSSPTree,
                           const unsigned i) noexcept
    requires (numThreads > 1)
    {
        DynamicGraphType<GraphType> tmpGraph(graph, result[i], reverseSSSPTree);

#pragma omp parallel for num_threads(numThreads), schedule(dynamic, 1), default(none), shared(t, i, reverseSSSPTree), firstprivate(tmpGraph)
        for(unsigned deviationNodeIndex = result[i].deviationNodeIndex; deviationNodeIndex < result[i].path.size() - 1; deviationNodeIndex++)
        {
            tmpGraph.advanceDeviationNodeIndex(deviationNodeIndex);

            if constexpr(attemptSSSPSkip > 0)
            {
                if(makeSSSPSkipAttempt<SSSPTreeType, attemptSSSPSkip >= 2, skipByLength>(i, deviationNodeIndex, tmpGraph, *reverseSSSPTree))
                    continue;
            }
            else
            {
                // if skips are not attempted by guessing the shortest deviation at least check if the deviation node
                // has any out edges left and skip SSSP construction and initialisation if so.
                auto outEdges = tmpGraph.template getOutEdges<true>(result[i].path[deviationNodeIndex]);

                // no valid out edge exists => go to the next deviation node
                if(outEdges.begin() == outEdges.end())
                    continue;
            }

            computeDeviationSuffix(tmpGraph, i, deviationNodeIndex, t);
        }
    }

    void computeDeviationSuffix(const DynamicGraphType<GraphType>& tmpGraph, const unsigned i,
                                const unsigned deviationNodeIndex, const NodeType t) noexcept
    {
        SSSPAlgorithmType<DynamicGraphType<GraphType>> tmpSSSP(tmpGraph);
        if constexpr(SSSPAlgorithmType<DynamicGraphType<GraphType>>::earlyStopping)
        {
            const auto prefixLength = deviationNodeIndex == 0
                                      ? static_cast<DistanceType>(0)
                                      : graph.getPathLength(std::span(result[i].path.begin(), deviationNodeIndex + 1));
            tmpSSSP.compute(result[i].path[deviationNodeIndex], t, candidates.getPathLengthLimit() - prefixLength);
        }
        else
        {
            tmpSSSP.compute(result[i].path[deviationNodeIndex], t, Distance<DistanceType>::max);
        }

        if constexpr(collectStatistics)
        {
            stats.storeNumExploredNodes(i, tmpSSSP.getSSSPTree()->getNumDiscoveredNodes());

            if constexpr(usesYellowGraph)
                stats.storeNumExpressEdges(i, tmpGraph.getNumExpressEdges());
        }

        if(!tmpSSSP.pathFound())
        {
            // we can assume that a path is only not found if it gets too long.
            stats.countSsspStoppedEarly();
            return;
        }

        auto pathPrefix = std::span(result[i].path.begin(), deviationNodeIndex + 1);
        auto pathSuffix = tmpSSSP.getPath();

        fixExpressEdges(pathSuffix, tmpGraph);

        addCandidatePath(i, deviationNodeIndex, pathPrefix + pathSuffix);
    }

    /**
     * Returns true if the SSSP computation can be skipped and false otherwise.
     */
    template<template<typename> typename SSSPTreeType, bool useSecondShortestPath, bool skipByLength>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    bool makeSSSPSkipAttempt(const unsigned parentPathId, const unsigned deviationNodeIndex,
                             const DynamicGraphType<GraphType>& tmpGraph,
                             const SSSPTreeType<DistanceType>& reverseSSSPTree) noexcept
    {
        auto outEdges = tmpGraph.template getOutEdges<true>(result[parentPathId].path[deviationNodeIndex]);
        auto lightestOutEdge = outEdges.begin();

        // no valid out edge exists => go to the next deviation node
        if(lightestOutEdge == outEdges.end())
            return true;

        if constexpr(usesYellowGraph)
        {
            // yellow graph uses express edges, so the guessed candidate can only be valid
            // if the target node of the edge is the actual target.
            if(lightestOutEdge->getTarget() == result[parentPathId].path.back())
            {
                const PathType candidate = std::span(result[parentPathId].path.begin(), deviationNodeIndex + 1)
                                           + tmpGraph.expressEdgeToPath(result[parentPathId].path[deviationNodeIndex]);

                assert(isPathSimple(candidate));

                addCandidatePath(parentPathId, deviationNodeIndex, candidate);
                stats.countSSCShortestDeviation();
                return true;
            }
            else    // if not, the target of the lightest edge is yellow and thus no express edge is generated
            {
                // do no additional work here if none of the both optimization are on.
                if constexpr(!skipByLength && !useSecondShortestPath)
                    return false;

                const PathType candidate = std::span(result[parentPathId].path.begin(), deviationNodeIndex + 1)
                                         + reverseSSSPTree.getReversePath(lightestOutEdge->getTarget());

                if constexpr(skipByLength)
                {
                    if(candidates.getPathLengthLimit() < Distance<DistanceType>::max &&
                       candidates.getPathLengthLimit() < graph.getPathLength(candidate))
                    {
                        stats.countSSCShortestDeviationLength();
                        return true;
                    }
                }

                if constexpr(useSecondShortestPath)
                {
                    return skipViaSecondShortestPath<skipByLength, SSSPTreeType>(tmpGraph, reverseSSSPTree, candidate, deviationNodeIndex, parentPathId);
                }
            }
        }
        else
        {
            NodeType start = lightestOutEdge->getTarget();
            DistanceType minDistance = lightestOutEdge->getWeight() + (guided ? 0 : reverseSSSPTree.getDistance(start));

            // If the graph is NOT guided, the smallest sum of edge and shortest path needs to be found first.
            if constexpr(!guided)
            {
                for(++lightestOutEdge; lightestOutEdge != outEdges.end(); ++lightestOutEdge)
                {
                    if(minDistance > lightestOutEdge->getWeight() + reverseSSSPTree.getDistance(lightestOutEdge->getTarget()))
                    {
                        minDistance = lightestOutEdge->getWeight() + reverseSSSPTree.getDistance(lightestOutEdge->getTarget());
                        start = lightestOutEdge->getTarget();
                    }
                }
            }

            const PathType candidate = std::span(result[parentPathId].path.begin(), deviationNodeIndex + 1)
                                       + reverseSSSPTree.getReversePath(start);

            if(isPathSimple(candidate))
            {
                addCandidatePath(parentPathId, deviationNodeIndex, candidate);
                stats.countSSCShortestDeviation();
                return true;
            }
            else
            {
                if constexpr(skipByLength)
                {
                    if(candidates.getPathLengthLimit() < Distance<DistanceType>::max &&
                       candidates.getPathLengthLimit() < graph.getPathLength(candidate))
                    {
                        stats.countSSCShortestDeviationLength();
                        return true;
                    }
                }

                if constexpr(useSecondShortestPath)
                {
                    return skipViaSecondShortestPath<skipByLength, SSSPTreeType>(tmpGraph, reverseSSSPTree, candidate, deviationNodeIndex, parentPathId);
                }
            }
        }

        return false;
    }

    template<bool skipByLength, template<typename> typename SSSPTreeType>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    bool skipViaSecondShortestPath(const DynamicGraphType<GraphType>& tmpGraph,
                                   const SSSPTreeType<DistanceType>& reverseSSSPTree,
                                   const PathType& firstShortestPath,
                                   const unsigned deviationNodeIndex,
                                   const unsigned parentPathId) noexcept
    {
        // at this point, firstShortestPath cannot be loopless.
        const auto firstIndexOfLoop = getFirstNonUniqueNode(firstShortestPath, deviationNodeIndex);
        assert(firstIndexOfLoop > 0);

        auto prefixLength = graph.getPathLength(std::span(firstShortestPath.begin(), deviationNodeIndex + 1));

        auto secondDeviationNodeIndex = deviationNodeIndex;
        auto [deviationNeighbor, secondShortestPathLength] = findSecondShortestPathNeighbor<SSSPTreeType, true>(
            tmpGraph, reverseSSSPTree, firstShortestPath[deviationNodeIndex], firstShortestPath[deviationNodeIndex + 1]
        );

        secondShortestPathLength += prefixLength;

        for(unsigned i = deviationNodeIndex + 1; i < firstIndexOfLoop; i++)
        {
            prefixLength += graph.getEdgeWeight(firstShortestPath[i - 1], firstShortestPath[i]);
            // find shortest out-edge that is not used in the first shortest path
            auto [newDeviationNeighbor, length] = findSecondShortestPathNeighbor<SSSPTreeType, false>(
                tmpGraph, reverseSSSPTree, firstShortestPath[i], firstShortestPath[i + 1]
            );

            length += prefixLength;

            if(length < secondShortestPathLength)
            {
                secondShortestPathLength = length;
                deviationNeighbor = newDeviationNeighbor;
                secondDeviationNodeIndex = i;
            }
        }

        if(secondShortestPathLength < Distance<DistanceType>::max)
        {
            const PathType candidate = std::span(firstShortestPath.begin(), secondDeviationNodeIndex + 1)
                                       + reverseSSSPTree.getReversePath(deviationNeighbor);

            if(isPathSimple(candidate))
            {
                stats.countSSCSecondShortestDeviation();
                addCandidatePath(parentPathId, deviationNodeIndex, candidate);
                return true;
            }
            else
            {
                if constexpr(skipByLength)
                {
                    if(candidates.getPathLengthLimit() < Distance<DistanceType>::max &&
                       candidates.getPathLengthLimit() < secondShortestPathLength)
                    {
                        stats.countSSCSecondShortestDeviationLength();
                        return true;
                    }
                }
            }
        }

        return false;
    }

    template<template<typename> typename SSSPTreeType, bool fromDeviationNode>
    requires SSSPTreeConcept<SSSPTreeType<DistanceType>>
    std::tuple<NodeType, DistanceType> findSecondShortestPathNeighbor(const DynamicGraphType<GraphType>& tmpGraph,
                                       [[maybe_unused]] const SSSPTreeType<DistanceType>& reverseSSSPTree,
                                       const NodeType source, const NodeType forbiddenNeighbor) const noexcept
    {
        // In case of a yellow graph, this turns off express edges
        // and allows to use the yellow graph as a normal guided graph
        constexpr bool noExpressEdges = true;
        auto outEdges = tmpGraph.template getOutEdges<fromDeviationNode, noExpressEdges>(source);

        if constexpr(guided)
        {
            auto lightestOutEdge = outEdges.begin();

            while(lightestOutEdge != outEdges.end() && lightestOutEdge->getTarget() == forbiddenNeighbor)
            {
                ++lightestOutEdge;
            }

            if(lightestOutEdge == outEdges.end())
                return std::tie(nullNode, Distance<DistanceType>::max);

            return std::make_tuple(lightestOutEdge->getTarget(), lightestOutEdge->getWeight());
        }
        else
        {
            DistanceType minDistance = Distance<DistanceType>::max;
            NodeType closestNeighbor = nullNode;

            for(const auto& edge : outEdges)
            {
                if(edge.getTarget() == forbiddenNeighbor)
                {
                    continue;
                }

                auto distance = edge.getWeight() + reverseSSSPTree.getDistance(edge.getTarget());
                if(minDistance > distance)
                {
                    minDistance = distance;
                    closestNeighbor = edge.getTarget();
                }
            }

            return std::tie(closestNeighbor, minDistance);
        }
    }

    void fixExpressEdges(PathType& pathSuffix, const DynamicGraphType<GraphType>& tmpGraph) const noexcept
    {
        if constexpr(usesYellowGraph)
        {
            // the suffix ends in node t, the node before t is the node that can needs to be
            // matched to the express edge/path.
            auto expressPath = tmpGraph.expressEdgeToPath(pathSuffix[pathSuffix.size() - 2]);
            // merge the express path into the path suffix.
            pathSuffix = std::span(pathSuffix.begin(), pathSuffix.size() - 1) + expressPath;
        }
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

    /**
     * Returns the index of the first node that introduces a loop, zero if no loop is detected.
     */
    [[nodiscard]] unsigned getFirstNonUniqueNode(const PathType& path, const unsigned deviationIndex) const noexcept
    {
        const auto prefix = std::span(path.begin(), deviationIndex + 1);
        for(unsigned i = deviationIndex + 1; i < path.size(); i++)
        {
            if(std::ranges::find(prefix, path[i]) != prefix.end())
                return i;
        }

        return 0;
    }

public:
    [[nodiscard]] const std::vector<KSPPath<DistanceType>>& getResults() const noexcept
    {
        return result;
    }

    [[nodiscard]] const Statistics<true>& getStats() const noexcept requires collectStatistics
    {
        return stats;
    }

    void printStatisticalData(const unsigned attemptSSSPSkip, const bool skipByLength, const size_t seed, const unsigned run) const noexcept
    {
        JSON json;
        json.add("algo", static_cast<std::string>(usesYellowGraph ? "feng" : "yen"));
        json.add("earlyStopping", SSSPAlgorithmType<DynamicGraphType<GraphType>>::earlyStopping);
        json.add("guided", guided);
        json.add("attemptSSSPSkips", attemptSSSPSkip);
        json.add("skipByLength", skipByLength);
        if constexpr(DeltaSteppingGraphConcept<GraphType>)
            json.add("delta", graph.getDelta());
        json.add("numNodes", graph.getNumNodes());
        json.add("numEdges", graph.getNumEdges());
        json.add("seed", seed);
        json.add("run", run);
        json.add("paths", result);
        json.add("stats", stats.template getJSON<usesYellowGraph>());

        Print::data() << json;
    }
};

template<GraphConcept GraphType, template<typename> typename SSSPAlgorithmType, unsigned numThreads = 1, bool collectStatistics = false>
using FengsAlgorithm = KSPAlgorithm<GraphType, SSSPAlgorithmType, YellowGraph, true, numThreads, collectStatistics>;

template<GraphConcept GraphType, template<typename> typename SSSPAlgorithmType, bool guided, unsigned numThreads = 1, bool collectStatistics = false>
using YensAlgorithm = KSPAlgorithm<GraphType, SSSPAlgorithmType, DynamicGraph, guided, numThreads, collectStatistics>;

#endif //KSP_YENSALGORITHM_H
