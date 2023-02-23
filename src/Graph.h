#ifndef KSP_GRAPH_H
#define KSP_GRAPH_H

#include "concepts/graph.h"
#include "concepts/sssp.h"

#include <algorithm>
#include <cassert>
#include <ranges>
#include <span>
#include <type_traits>
#include <vector>


class BasicEdge
{
    static_assert(std::is_unsigned<NodeType>::value, "NodeType needs to be unsigned.");

private:
    NodeType target;

public:
    explicit BasicEdge(const NodeType target) noexcept : target(target)
    {}

    [[nodiscard]] NodeType getTarget() const noexcept
    {
        return target;
    }

    void setTarget(const NodeType newTarget) noexcept
    {
        target = newTarget;
    }
};

template<bool directed>
class DirectedOption;

template<>
class DirectedOption<true>
{
public:
    [[nodiscard]] EdgeIdType getBackwardEdgeId() const noexcept
    {
        assert(false);
        return nullEdge;
    }

    void setBackwardEdgeId(EdgeIdType) noexcept
    {
        assert(false);
    }
};

template<>
class DirectedOption<false>
{
public:
    [[nodiscard]] EdgeIdType getBackwardEdgeId() const noexcept
    {
        return backwardEdgeId;
    }

    void setBackwardEdgeId(EdgeIdType newBackwardEdgeId) noexcept
    {
        backwardEdgeId = newBackwardEdgeId;
    }

private:
    EdgeIdType backwardEdgeId = nullEdge;
};

template<typename T>
class WeightOption
{
public:
    static_assert(WeightConcept<T>, "Invalid weight type");
    using WeightType = T;

private:
    WeightType weight;

public:
    explicit WeightOption(WeightType weight) noexcept : weight(weight)
    {}

    [[nodiscard]] WeightType getWeight() const noexcept
    {
        return weight;
    }

    void setWeight(WeightType newWeight) noexcept
    {
        weight = newWeight;
    }
};

template<>
class WeightOption<void>
{
public:
    using WeightType = float;

    [[nodiscard]] constexpr WeightType getWeight() const noexcept // NOLINT(readability-convert-member-functions-to-static) need to match the signature of the weighted version.
    {
        return 1.0f;
    }
};

/**
 *
 * @tparam WT "void" is interpreted as unweighted, numeric types are used as Weight type, all other types are invalid.
 * @tparam directed
 */
template<typename WT, bool directed>
class Edge
    : public BasicEdge
    , public WeightOption<WT>
    , public DirectedOption<directed>
{
public:
    static constexpr bool isWeighted = !std::same_as<WT, void>;
    static constexpr bool isDirected = directed;

    using WeightType = typename WeightOption<WT>::WeightType;

    explicit Edge(NodeType target) noexcept requires (!isWeighted)
      : BasicEdge(target)
    {}

    explicit Edge(NodeType target, WeightType weight) noexcept requires isWeighted
      : BasicEdge(target), WeightOption<WT>(weight)
    {}

    auto operator<=>(const Edge<WT, directed>& rhs) const noexcept
    {
        if constexpr(isWeighted)
            return this->getWeight() <=> rhs.getWeight();
        else
            return 0;
    }
};

/**
 * This is the generic data structure for a graph.
 * Since it is designed to work for k shortest path algorithms it supports functions
 * to temporarily remove nodes and edges if the edges support it.
 *
 * @tparam EdgeType
 * @tparam readyForDeltaStepping
 */
template<EdgeConcept ET, bool readyForDeltaStepping = true>
class Graph
{
public:
    using EdgeType = ET;

    static constexpr bool isDirected = EdgeType::isDirected;
    static constexpr bool isWeighted = EdgeType::isWeighted;

    using EdgeListType = std::vector<EdgeType>;
    using EdgeIdType = typename EdgeListType::size_type;
    using EdgeIdListType = std::vector<EdgeIdType>;

    using WeightType = typename EdgeType::WeightType;
    using DistanceType = WeightType;

    const NodeType _numNodes;
    const EdgeIdListType _outEdgesBeginnings;
    const EdgeListType _edges;
    const WeightType _delta = 1.0;
    const WeightType _heaviestWeight = 1.0;

public:

    /**
     * Constructor used for Graphs not used for delta stepping.
     * @param numNodes
     * @param outEdgesBegin Edge Id where neighbors begin. Light neighbors of node u start at index 2u, heavy neighbors
     * at 2u+1. The weight of light edges is at most delta. The list ends with a sentinel.
     * @param edges Edges (u,v) grouped by tail nodes u and sorted by weight.
     */
    Graph(const NodeType numNodes, EdgeIdListType& outEdgesBeginnings, EdgeListType& edges, const WeightType delta, const WeightType heaviestWeight) noexcept requires readyForDeltaStepping
      : _numNodes(numNodes), _outEdgesBeginnings(std::move(outEdgesBeginnings)), _edges(std::move(edges)), _delta(delta), _heaviestWeight(heaviestWeight)
    {
        assert(_outEdgesBeginnings.size() == 2 * _numNodes + 1);
        assert(_outEdgesBeginnings[getOutEdgeBeginningId(_numNodes)] == _edges.size());
        assert(heaviestWeight >= delta);
    }

    /**
     * Constructor used for Graphs not used for delta stepping.
     * @param numNodes
     * @param outEdgesBegin Edge Ids where neighbors of a node begin. Ends with a sentinel.
     * @param edges Edges (u,v) grouped by tail nodes u.
     */
    Graph(const NodeType numNodes, EdgeIdListType& outEdgesBeginnings, EdgeListType& edges) noexcept requires (!readyForDeltaStepping)
      : _numNodes(numNodes), _outEdgesBeginnings(std::move(outEdgesBeginnings)), _edges(std::move(edges))
    {
        assert(_outEdgesBeginnings.size() == _numNodes + 1);
        assert(_outEdgesBeginnings[getOutEdgeBeginningId(_numNodes)] == _edges.size());
    }

    static Graph<EdgeType, readyForDeltaStepping> create(const std::vector<std::vector<Edge<float, true>>>& adjacencyList) noexcept requires (!readyForDeltaStepping)
    {
        const auto numNodes = static_cast<NodeType>(adjacencyList.size());
        std::vector<EdgeType> outEdges;
        EdgeIdListType outEdgesBeginnings(numNodes + 1, 0);

        for(NodeType i = 0; i < numNodes; i++)
        {
            outEdgesBeginnings[i + 1] = outEdgesBeginnings[i] + adjacencyList[i].size();
            for(const auto& edge : adjacencyList[i])
            {
                if constexpr(isWeighted)
                    outEdges.emplace_back(edge.getTarget(), static_cast<WeightType>(edge.getWeight()));
                else
                    outEdges.template emplace_back(edge.getTarget());
            }
        }

        assert(outEdges.size() == outEdgesBeginnings[numNodes]);

        if(!isDirected)
        {
            // todo link backward edges
            assert(false);
        }

        return {numNodes, outEdgesBeginnings, outEdges};
    }

public:
    [[nodiscard]] Graph<EdgeType, readyForDeltaStepping> computeReverseGraph() const noexcept
    {
        return _computeReverseGraph<false>();
    }

    [[nodiscard]] Graph<EdgeType, readyForDeltaStepping> computeReverseGraph(const NodeType s, const NodeType t) const noexcept
    {
        assert(s != nullNode);
        assert(t != nullNode);
        return _computeReverseGraph<true>(s, t);
    }

private:
    template<bool removeInOutEdges = false>
    [[nodiscard]] Graph<EdgeType, readyForDeltaStepping> _computeReverseGraph(const NodeType s = nullNode, const NodeType t = nullNode) const noexcept
    {
        if constexpr(!isDirected)
            return this;

        EdgeListType reverseEdges;
        reverseEdges.reserve(_edges.size());
        EdgeIdListType reverseOutEdgesBeginnings(_outEdgesBeginnings.size(), 0);   // includes sentinel already

        std::vector<std::pair<NodeType, EdgeType>> tmpEdgeList;
        tmpEdgeList.reserve(_edges.size());

        for(NodeType u = 0; u < _numNodes; u++)
        {
            if constexpr(removeInOutEdges)
            {
                if(u == t)  // skip all out-edges of the target node
                    continue;
            }

            for(const auto& e : getOutEdges(u))
            {
                if constexpr(removeInOutEdges)
                {
                    if(e.getTarget() == s)  // skip all in-edges of the source node
                        continue;
                }

                if constexpr(isWeighted)
                    tmpEdgeList.emplace_back(e.getTarget(), EdgeType(u, e.getWeight()));
                else
                    tmpEdgeList.emplace_back(e.getTarget(), EdgeType(u));
            }
        }

        std::ranges::sort(tmpEdgeList, [](const auto& lhs, const auto& rhs){ return lhs.first < rhs.first; });

        NodeType currentNode = 0;

        const auto tmpDelta = _delta;

        for(const auto& [u, e] : tmpEdgeList)
        {
            if(u != currentNode)
            {
                assert(u > currentNode);
                if constexpr(!removeInOutEdges)
                    assert(u == currentNode + 1);

                sortOutEdges(currentNode, tmpDelta, reverseEdges, reverseOutEdgesBeginnings);

                // make sure to set reverseOutEdgesBeginnings for nodes without out-edges too
                for(NodeType j = currentNode + 1; j <= u; j++)
                {
                    reverseOutEdgesBeginnings[getOutEdgeBeginningId(j)] = reverseEdges.size();
                    if constexpr(readyForDeltaStepping)
                        reverseOutEdgesBeginnings[getHeavyOutEdgeBeginningId(j)] = reverseEdges.size();
                }

                currentNode = u;
            }

            reverseEdges.push_back(e);
        }

        // set IDs of out edge beginnings of nodes the last nodes that have no out edges
        // and thus were not hit by the loop. This includes the sentinel.
        for(size_t i = getOutEdgeBeginningId(currentNode + 1); i < reverseOutEdgesBeginnings.size(); i++)
        {
            reverseOutEdgesBeginnings[i] = reverseEdges.size();
        }

        // sort the edges of the final node.
        sortOutEdges(currentNode, tmpDelta, reverseEdges, reverseOutEdgesBeginnings);

        if constexpr(!removeInOutEdges)
            assert(reverseEdges.size() == _edges.size());

        if constexpr(readyForDeltaStepping)
            return Graph<EdgeType, true>(_numNodes, reverseOutEdgesBeginnings, reverseEdges, tmpDelta, _heaviestWeight);
        else
            return Graph<EdgeType, false>(_numNodes, reverseOutEdgesBeginnings, reverseEdges);
    }

    static void sortOutEdges(const NodeType node, const WeightType tmpDelta, EdgeListType& reverseEdges, EdgeIdListType& reverseOutEdgesBeginnings) noexcept
    {
        const auto edgesBegin = reverseEdges.begin() + static_cast<EdgeIdDifferenceType>(reverseOutEdgesBeginnings[getOutEdgeBeginningId(node)]);
        const auto edgesEnd = reverseEdges.end();
        std::sort(edgesBegin, edgesEnd);

        if constexpr(readyForDeltaStepping)
        {
            const auto tmpHeavyEdgesBegin = std::partition(edgesBegin, edgesEnd, [tmpDelta](const auto& edge){ return edge.getWeight() < tmpDelta; });
            assert(edgesBegin <= tmpHeavyEdgesBegin);
            assert(tmpHeavyEdgesBegin <= edgesEnd);
            reverseOutEdgesBeginnings[getHeavyOutEdgeBeginningId(node)] = static_cast<EdgeIdType>(std::distance(reverseEdges.begin(), tmpHeavyEdgesBegin));
        }
    }

public:
    template<template<typename> typename SSSPAlgorithmType>
    [[nodiscard]] auto precomputeSSSPGuiding(const NodeType source, const NodeType target) const noexcept
    {
        return precomputeSSSPGuiding<SSSPAlgorithmType>(source, target, _delta);
    }

    /*
     * Compute the reverse sssp tree rooted at t and recalculate edge weights to guide SSSP algorithms
     * into the right direction
     */
    template<template<typename> typename SSSPAlgorithmType>
    requires SSSPAlgoConcept<SSSPAlgorithmType<Graph<ET, readyForDeltaStepping>>>
    [[nodiscard]] auto precomputeSSSPGuiding(const NodeType source, const NodeType target, const WeightType newDelta) const noexcept
    {
        using NewWeightType = std::conditional_t<isWeighted, WeightType, std::uint32_t>;
        using NewEdgeType = Edge<NewWeightType, isDirected>;
        auto reverseGraph = computeReverseGraph(source, target);

        SSSPAlgorithmType ssspAlgo(reverseGraph);
        ssspAlgo.compute(target);
        auto ssspTree = ssspAlgo.getSSSPTree();

        // update edge weights
        std::vector<EdgeIdType> outEdgesBeginnings(_outEdgesBeginnings.size(), 0);
        std::vector<NewEdgeType> edges;
        edges.reserve(_edges.size());

        NewWeightType newHeaviestWeight = 0;

        for(NodeType node = 0; node < _numNodes; node++)
        {
            const auto nodeDistance = ssspTree->getDistance(node);

            // skip all out-edges if node is not connected to the target or is the target itself
            if(nodeDistance != Distance<DistanceType>::max && node != target)
            {
                for(const auto& edge : getOutEdges(node))
                {
                    const auto edgeTarget = edge.getTarget();
                    const auto targetDistance = ssspTree->getDistance(edgeTarget);

                    // skip all in-edges of s and edges that cannot lead to the target
                    if(targetDistance == Distance<DistanceType>::max || edgeTarget == source)
                        continue;

                    const NewWeightType newWeight = edge.getWeight() + targetDistance - nodeDistance;

                    assert(newWeight >= 0);

                    newHeaviestWeight = std::max(newHeaviestWeight, newWeight);
                    edges.emplace_back(edgeTarget, newWeight);
                }
            }

            outEdgesBeginnings[getOutEdgeBeginningId(node + 1)] = edges.size();

            const auto outEdgesBegin = edges.begin() + static_cast<EdgeIdDifferenceType>(outEdgesBeginnings[getOutEdgeBeginningId(node)]);
            const auto outEdgesEnd = edges.end();

            if constexpr(readyForDeltaStepping)
            {
                if(outEdgesBegin == outEdgesEnd)
                {
                    outEdgesBeginnings[getHeavyOutEdgeBeginningId(node)] = outEdgesBeginnings[getOutEdgeBeginningId(node + 1)];
                }
                else
                {
                    const auto firstHeavyEdge = std::partition(outEdgesBegin, outEdgesEnd,
                                                               [newDelta](const auto& edge) {
                                                                   return edge.getWeight() < newDelta;
                                                               });
                    outEdgesBeginnings[getHeavyOutEdgeBeginningId(node)] = static_cast<EdgeIdType>(std::distance(edges.begin(), firstHeavyEdge));
                }

                assert(outEdgesBeginnings[getOutEdgeBeginningId(node)] <= outEdgesBeginnings[getHeavyOutEdgeBeginningId(node)]);
            }

            std::sort(outEdgesBegin, outEdgesEnd);
        }

        if constexpr(readyForDeltaStepping)
            return std::tuple(Graph<NewEdgeType, readyForDeltaStepping>(_numNodes, outEdgesBeginnings, edges, newDelta, newHeaviestWeight), ssspTree);
        else
            return std::tuple(Graph<NewEdgeType, readyForDeltaStepping>(_numNodes, outEdgesBeginnings, edges), ssspTree);
    }

    [[nodiscard]] NodeType getNumNodes() const noexcept
    {
        return _numNodes;
    }

    [[nodiscard]] EdgeIdType getNumEdges() const noexcept
    {
        return _edges.size();
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] std::span<const EdgeType> getOutEdges(const NodeType node) const noexcept
    {
        return std::span(_edges.data() + getNeighborBeginId(node), getNumNeighbors(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getOutEdgesBegin(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getNeighborBeginId(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getOutEdgesEnd(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getNeighborEndId(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] std::span<const EdgeType> getLightOutEdges(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return std::span(_edges.data() + getLightNeighborBeginId(node), getNumLightNeighbors(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getLightOutEdgesBegin(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getLightNeighborBeginId(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getLightOutEdgesEnd(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getLightNeighborEndId(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] std::span<const EdgeType> getHeavyOutEdges(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return std::span(_edges.data() + getHeavyNeighborBeginId(node), getNumHeavyNeighbors(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getHeavyOutEdgesBegin(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getHeavyNeighborBeginId(node));
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] auto getHeavyOutEdgesEnd(const NodeType node) const noexcept
    {
        return _edges.cbegin() + static_cast<typename std::vector<EdgeType>::const_iterator::difference_type>(getHeavyNeighborEndId(node));
    }

    [[nodiscard]] WeightType getDelta() const noexcept requires readyForDeltaStepping
    {
        return _delta;
    }

    [[nodiscard]] WeightType getHeaviestWeight() const noexcept requires readyForDeltaStepping
    {
        return _heaviestWeight;
    }

    /**
     * This function assumes that each edge on the path exists.
     */
    [[nodiscard]] DistanceType getPathLength(const std::span<const NodeType> path) const noexcept
    {
        DistanceType length = 0;

        for(unsigned i = 0; i < path.size() - 1; i++)
            length += getEdgeWeight(path[i], path[i + 1]);

        return length;
    }

    /**
     * This function assumes that the edge (u,v) exists.
     */
    [[nodiscard]] WeightType getEdgeWeight(const NodeType u, const NodeType v) const noexcept
    {
        for(const auto& e : getOutEdges(u))
        {
            if(e.getTarget() == v)
                return e.getWeight();
        }

        assert(false);
        return Weight<WeightType>::max;
    }

private:
    [[nodiscard]] static constexpr NodeType getOutEdgeBeginningId(const NodeType node) noexcept
    {
        if constexpr(readyForDeltaStepping)
            return 2 * node;
        else
            return node;
    }

    [[nodiscard]] static constexpr NodeType getHeavyOutEdgeBeginningId(const NodeType node) noexcept requires readyForDeltaStepping
    {
        return getOutEdgeBeginningId(node) + 1;
    }

    [[nodiscard]] EdgeIdType getNeighborBeginId(const NodeType node) const noexcept
    {
        return _outEdgesBeginnings[getOutEdgeBeginningId(node)];
    }

    [[nodiscard]] EdgeIdType getLightNeighborBeginId(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return getNeighborBeginId(node);
    }

    [[nodiscard]] EdgeIdType getLightNeighborEndId(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return _outEdgesBeginnings[getHeavyOutEdgeBeginningId(node)];
    }

    [[nodiscard]] EdgeIdType getHeavyNeighborBeginId(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return getLightNeighborEndId(node);
    }

    [[nodiscard]] EdgeIdType getHeavyNeighborEndId(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return getNeighborEndId(node);
    }

    [[nodiscard]] EdgeIdType getNeighborEndId(const NodeType node) const noexcept
    {
        return getNeighborBeginId(node + 1);
    }

    [[nodiscard]] EdgeIdType getNumNeighbors(const NodeType node) const noexcept
    {
        return getNeighborEndId(node) - getNeighborBeginId(node);
    }

    [[nodiscard]] EdgeIdType getNumLightNeighbors(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return getLightNeighborEndId(node) - getLightNeighborBeginId(node);
    }

    [[nodiscard]] EdgeIdType getNumHeavyNeighbors(const NodeType node) const noexcept requires readyForDeltaStepping
    {
        return getHeavyNeighborEndId(node) - getHeavyNeighborBeginId(node);
    }

};


#endif //KSP_GRAPH_H
