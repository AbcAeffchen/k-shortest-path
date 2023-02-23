#ifndef KSP_GRAPH_CONCEPT_H
#define KSP_GRAPH_CONCEPT_H

#include <algorithm>
#include <cassert>
#include <concepts>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>
#include <vector>

/* ******************************************
 * Most basic node and edge definitions.
 ****************************************** */

using NodeType = uint32_t;
constexpr NodeType nullNode = std::numeric_limits<NodeType>::max(); // maximal node id
//constexpr NodeType nullNode = std::numeric_limits<NodeType>::max() >> 1; // maximal node id without using the most significant bit

using EdgeIdType = size_t;
using EdgeIdDifferenceType = std::make_signed<EdgeIdType>::type;
constexpr EdgeIdType nullEdge = std::numeric_limits<EdgeIdType>::max();

using PathType = std::vector<NodeType>;

/**
 * Join two paths. It is assumed that either
 * 1) The last node of the prefix is the same as the first node of the suffix path. Then the node is only included once.
 * 2) The last node of the prefix is different from the first node of the suffix path. In this case both nodes are
 *    included and it is assumed that the edge exists in the graph.
 */
template<typename T, typename S>
requires std::ranges::range<T> && std::ranges::range<S>
inline PathType operator+(const T& prefix, const S& suffix)
{
    unsigned offset = prefix.back() == suffix.front() ? 1 : 0;
    PathType result;
    result.reserve(prefix.size() + suffix.size() - offset);
    result.insert(result.end(), prefix.begin(), prefix.end());
    result.insert(result.end(), suffix.begin() + offset, suffix.end());

    return result;
}

inline bool isPathSimple(PathType path)
{
    std::ranges::sort(path);
    for(size_t i = 1; i < path.size(); i++)
    {
        if(path[i - 1] == path[i])
            return false;
    }

    return true;
}

/*
 * Prints a path in a readable form for debugging.
 */
inline std::ostream& operator<<(std::ostream& os, const PathType& path)
{
    os << path[0];
    for(const auto node : path | std::views::drop(1))
        os << " -> " << node;

    return os;
}

/* ******************************************
 * Graph and graph related concepts
 ****************************************** */

template<typename T>
concept RangeConcept = requires (T range) {
    range.begin();
    range.end();
};

template<typename T>
concept WeightConcept = std::is_arithmetic_v<T> && std::numeric_limits<T>::has_infinity && !std::same_as<T, bool>;

template<WeightConcept WeightType>
struct Weight
{
    static constexpr WeightType min = static_cast<WeightType>(0);
    static constexpr WeightType max = std::numeric_limits<WeightType>::infinity();
};

template<WeightConcept WeightType>
using Distance = Weight<WeightType>;

template<typename T>
concept EdgeConcept = requires (T edge, T edge2, EdgeIdType eid) {
    requires std::same_as<const bool, decltype(T::isWeighted)>;
    requires std::same_as<const bool, decltype(T::isDirected)>;
    requires WeightConcept<typename T::WeightType>;
    { edge.getTarget() } -> std::same_as<NodeType>;
    { edge.getWeight() } -> std::same_as<typename T::WeightType>;
    { edge.getBackwardEdgeId() } -> std::same_as<EdgeIdType>;
    { edge.setBackwardEdgeId(eid) } -> std::same_as<void>;
    edge < edge2;
};

template<typename T>
concept GraphConcept = requires (T graph, NodeType node, PathType path) {
    requires std::same_as<const bool, decltype(T::isWeighted)>;
    requires std::same_as<const bool, decltype(T::isDirected)>;
    requires EdgeConcept<typename T::EdgeType>;
    requires std::integral<typename T::EdgeIdType>;
    requires WeightConcept<typename T::WeightType>;
    requires WeightConcept<typename T::DistanceType>;
    { graph.getNumNodes() } -> std::same_as<NodeType>;
    { graph.getNumEdges() } -> std::same_as<EdgeIdType>;
    { graph.getPathLength(path) } -> std::same_as<typename T::DistanceType>;
    { graph.getEdgeWeight(node, node) } -> std::same_as<typename T::DistanceType>;
    requires RangeConcept<decltype(graph.template getOutEdges<true>(node))>;  // allows std::vector, std::span
    requires RangeConcept<decltype(graph.template getOutEdges<false>(node))>;  // allows std::vector, std::span
};

template<typename T>
concept DeltaSteppingGraphConcept = requires (T graph, NodeType node)
{
    requires GraphConcept<T>;
    requires RangeConcept<decltype(graph.template getLightOutEdges<true>(node))>;
    requires RangeConcept<decltype(graph.template getLightOutEdges<false>(node))>;
    requires RangeConcept<decltype(graph.template getHeavyOutEdges<true>(node))>;
    requires RangeConcept<decltype(graph.template getHeavyOutEdges<false>(node))>;
    { graph.getDelta() } -> std::same_as<typename T::WeightType>;
    { graph.getHeaviestWeight() } -> std::same_as<typename T::WeightType>;
};

#endif //KSP_GRAPH_CONCEPT_H
