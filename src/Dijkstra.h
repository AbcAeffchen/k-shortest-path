#ifndef KSP_DIJKSTRA_H
#define KSP_DIJKSTRA_H

#include <cmath>

#include "concepts/sssp.h"
#include "PriorityQueue.h"
#include "SsspTree.h"

/**
 *
 * @tparam GraphType
 * @tparam stopEarly
 */
template<GraphConcept GraphType, template<typename> typename SSSPTreeType = SSSPTree>
requires SSSPTreeConcept<SSSPTreeType<typename GraphType::WeightType>>
class Dijkstra
{
public:
    using DistanceType = typename GraphType::WeightType;

private:
    enum class COMPUTE_MODE
    {
        FULL, STOP_WHEN_TARGET_IS_FOUND, STOP_WHEN_TARGET_IS_TOO_FAR_AWAY
    };

    const GraphType& graph;

    std::shared_ptr<SSSPTreeType<DistanceType>> ssspTree;
    NodeType _target = nullNode;

    // priority queue of tuples containing distance, node and predecessor leading to this node with this distance.
    PriorityQueue<DistanceType> pq;

    bool computationFinished = false;

public:
    explicit Dijkstra(const GraphType& graph)
      : graph(graph)
      , ssspTree(std::make_shared<SSSPTreeType<DistanceType>>(graph.getNumNodes()))
      , pq(static_cast<size_t>(std::sqrt(graph.getNumNodes())))  // sqrt is just a guess here. maybe there is a better value.
    {}

    /**
     * Computation stops as soon as the target is found or if the target is further away than the maximal allowed distance.
     */
    void compute(const NodeType start, const NodeType target, DistanceType maxDistance)
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        if(maxDistance == Distance<DistanceType>::max)
            _compute<COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND>(start, target);
        else
            _compute<COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY>(start, target, maxDistance);
    }

    /**
     * Stops as soon as the target is found.
     */
    void compute(const NodeType start, const NodeType target)
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        _compute<COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND>(start, target);
    }

    /**
     * Computes the full SSSP tree.
     */
    void compute(const NodeType start)
    {
        _compute<COMPUTE_MODE::FULL>(start);
        ssspTree->setChildLists();
    }

private:
    /**
     * Computes the shortest path from start to target.
     * @param start
     * @param target
     */
    template<COMPUTE_MODE cm>
    void _compute(const NodeType start, const NodeType target = nullNode, const DistanceType maxDistance = Distance<DistanceType>::max)
    {
        ssspTree->set(start, 0, nullNode);

        if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
        {
            if(start == target)
            {
                computationFinished = true;
                return;
            }
        }

        for(const auto& edge : graph.template getOutEdges<true>(start))
        {
            pq.template push<false>(edge.getWeight(), edge.getTarget(), start);
        }

        pq.update();

        while(!pq.empty())
        {
            const auto [newDistance, currentNode, predecessor] = pq.top();
            pq.pop();

            if constexpr(cm == COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY)
            {
                if(newDistance > maxDistance)
                    return;
            }

            // check if node is already settled and continue if so.
            if(newDistance >= ssspTree->getDistance(currentNode))
                continue;

            ssspTree->set(currentNode, newDistance, predecessor);

            // if stopEarly is active, we try to stop before unnecessary relaxing of out-edges
            if constexpr (cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
            {
                if(currentNode == target)
                    break;
            }

            // relax out-edges
            for(const auto& edge : graph.template getOutEdges<false>(currentNode))
            {
                pq.push(newDistance + edge.getWeight(), edge.getTarget(), currentNode);
            }
        }

        computationFinished = true;
    }

public:
    [[nodiscard]] bool pathFound() const
    {
        assert(_target != nullNode);
        return pathFound(_target);
    }

    [[nodiscard]] bool pathFound(const NodeType target) const
    {
        assert(_target == nullNode || _target == target);
        return computationFinished && ssspTree->getDistance(target) != Distance<DistanceType>::max;
    }

    [[nodiscard]] DistanceType getDistance() const
    {
        assert(_target != nullNode);
        return getDistance(_target);
    }

    [[nodiscard]] DistanceType getDistance(const NodeType target) const
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getDistance(target);
    }

    [[nodiscard]] PathType getPath() const
    {
        assert(_target != nullNode);
        return getPath(_target);
    }

    [[nodiscard]] PathType getPath(const NodeType target) const
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getPath(target);
    }

    std::shared_ptr<const SSSPTreeType<DistanceType>> getSSSPTree() const
    {
        return ssspTree;
    }
};


#endif //KSP_DIJKSTRA_H
