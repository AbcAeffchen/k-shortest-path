#ifndef KSP_SSSP_H
#define KSP_SSSP_H

#include "graph.h"

template<typename T>
concept SSSPTreeConcept = requires(T ssspTree, NodeType node, typename T::DistanceType distance) {
    WeightConcept<typename T::DistanceType>;
    { T(node) };
    { ssspTree.set(node, distance) } -> std::same_as<void>;
    { ssspTree.set(node, distance, node) } -> std::same_as<void>;
    { ssspTree.getDistance(node) } -> std::same_as<typename T::DistanceType>;
    { ssspTree.getPath(node) } -> std::same_as<PathType>;
};

template<typename T>
concept SSSPAlgoConcept = requires (T ssspAlgo, NodeType s, NodeType t, typename T::DistanceType weight) {
    requires WeightConcept<typename T::DistanceType>;
    { ssspAlgo.compute(s) } -> std::same_as<void>;
    { ssspAlgo.compute(s, t) } -> std::same_as<void>;
    { ssspAlgo.compute(s, t, weight) } -> std::same_as<void>;
    { ssspAlgo.pathFound() } -> std::same_as<bool>;
    { ssspAlgo.pathFound(t) } -> std::same_as<bool>;
    { ssspAlgo.getPath() } -> std::same_as<PathType>;
    { ssspAlgo.getPath(t) } -> std::same_as<PathType>;
    { ssspAlgo.getDistance() } -> std::same_as<typename T::DistanceType>;
    { ssspAlgo.getDistance(t) } -> std::same_as<typename T::DistanceType>;
    { ssspAlgo.getSSSPTree() };
};

#endif //KSP_SSSP_H
