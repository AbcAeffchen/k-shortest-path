#ifndef SRC_DYNAMICGRAPH
#define SRC_DYNAMICGRAPH

#include "concepts/graph.h"
#include "SsspTree.h"
#include "KspBasics.h"

/**
 * Lightweight wrapper to make a graph dynamic.
 * This wrapper assumes that only very few nodes and edges get removed.
 */
template<GraphConcept GraphType>
class DynamicGraph
{
public:
    static constexpr bool isDirected = GraphType::isDirected;
    static constexpr bool isWeighted = GraphType::isWeighted;
    static constexpr bool isDynamic = GraphType::isDynamic;

    using EdgeType = typename GraphType::EdgeType;
    using EdgeIdType = typename GraphType::EdgeIdType;
    using WeightType = typename GraphType::WeightType;
    using DistanceType = typename GraphType::DistanceType;

private:
    using ParentPathType = KSPPath<WeightType>;

public:
    template<bool fromSourceNode>
    class EdgeRange
    {
        using BaseEdgeListType = std::vector<EdgeType>;
        using BaseEdgeListIteratorType = typename BaseEdgeListType::const_iterator;

        class AdditionalIteratorFields
        {
        protected:
            const NodeType deviationNodeIndex;
            const ParentPathType& parentPath;

        public:
            AdditionalIteratorFields(NodeType deviationNodeIndex, const ParentPathType& parentPath) noexcept
                : deviationNodeIndex(deviationNodeIndex), parentPath(parentPath)
            {}
        };

        class NoAdditionalIteratorFields
        {};

    public:
        class Iterator : public std::conditional_t<fromSourceNode, AdditionalIteratorFields, NoAdditionalIteratorFields>
        {
            BaseEdgeListIteratorType currentEdge;
            const BaseEdgeListIteratorType end;
            const std::vector<bool>& removedNodes;

        public:
            Iterator(const BaseEdgeListIteratorType begin, const BaseEdgeListIteratorType end,
                     const NodeType deviationNodeIndex, const std::vector<bool>& removedNodes,
                     const ParentPathType& parentPath) noexcept requires fromSourceNode
                : AdditionalIteratorFields(deviationNodeIndex, parentPath)
                , currentEdge(begin), end(end)
                , removedNodes(removedNodes)
            {
                while(currentEdge != end && !currentEdgeIsInRange())
                    ++currentEdge;
            }

            Iterator(const BaseEdgeListIteratorType& begin, const BaseEdgeListIteratorType& end,
                     const std::vector<bool>& removedNodes) noexcept requires (!fromSourceNode)
                : currentEdge(begin), end(end)
                , removedNodes(removedNodes)
            {
                while(currentEdge != end && !currentEdgeIsInRange())
                    ++currentEdge;
            }

            inline void operator++() noexcept
            {
                do
                {
                    ++currentEdge;
                }
                while(currentEdge != end && !currentEdgeIsInRange());
            }

            inline const EdgeType& operator*() const noexcept
            {
                return *currentEdge;
            }

            inline auto operator->() const noexcept
            {
                return currentEdge;
            }

            inline friend bool operator!=(const Iterator& iterator, const BaseEdgeListIteratorType rhs) noexcept
            {
                return iterator.currentEdge != rhs;
            }

            inline friend bool operator==(const Iterator& iterator, const BaseEdgeListIteratorType rhs) noexcept
            {
                return iterator.currentEdge == rhs;
            }

        private:
            [[nodiscard]] bool currentEdgeIsInRange() const noexcept
            {
                if constexpr(fromSourceNode)
                {
                    if(currentEdge->getTarget() == this->parentPath.path[this->deviationNodeIndex + 1]
                       || (this->deviationNodeIndex == this->parentPath.deviationNodeIndex
                           && std::ranges::find(this->parentPath.forbiddenEdges, currentEdge->getTarget()) != this->parentPath.forbiddenEdges.end()))
                    {
                        return false;
                    }
                }

                return !removedNodes[currentEdge->getTarget()];
            }
        };

        const BaseEdgeListIteratorType allEdgesBegin;
        const BaseEdgeListIteratorType allEdgesEnd;
        const NodeType deviationNodeIndex;
        const std::vector<bool>& removedNodes;
        const ParentPathType& parentPath;

    public:
        EdgeRange(const NodeType deviationNodeIndex,
                  const BaseEdgeListIteratorType allEdgesBegin, const BaseEdgeListIteratorType allEdgesEnd,
                  const std::vector<bool>& removedNodes, const ParentPathType& parentPath) noexcept
            : allEdgesBegin(allEdgesBegin), allEdgesEnd(allEdgesEnd)
            , deviationNodeIndex(deviationNodeIndex)
            , removedNodes(removedNodes)
            , parentPath(parentPath)
        {}

        auto begin() const noexcept
        {
            if constexpr(fromSourceNode)
                return Iterator(allEdgesBegin, allEdgesEnd, deviationNodeIndex, removedNodes, parentPath);
            else
                return Iterator(allEdgesBegin, allEdgesEnd, removedNodes);
        }

        auto end() const noexcept
        {
            return allEdgesEnd;
        }
    };

private:
    const GraphType& graph;
    const ParentPathType& parentPath;
    NodeType deviationNodeIndex = 0;
    std::vector<bool> removedNodes;

public:
    DynamicGraph(const GraphType& graph, const KSPPath<typename GraphType::WeightType>& parentPath,
                 void const * const) noexcept
        : graph(graph), parentPath(parentPath), removedNodes(graph.getNumNodes(), false)
    {}

    DynamicGraph(const GraphType& graph, const KSPPath<typename GraphType::WeightType>& parentPath,
                 unsigned deviationNodeIndex, void const * const) noexcept
        : graph(graph), parentPath(parentPath), removedNodes(graph.getNumNodes(), false)
    {
        advanceDeviationNodeIndex(deviationNodeIndex);
    }

    // copy constructor needed for multi-threaded KSP.
    DynamicGraph(const DynamicGraph<GraphType>& dynamicGraph) = default;

    /*
     * Use this to recycle the object. For some reason it seems to be faster to create the object from scratch...
     */
    void advanceDeviationNodeIndex(const NodeType newDeviationNodeIndex) noexcept
    {
        assert(parentPath.deviationNodeIndex <= newDeviationNodeIndex && parentPath.path.size() - 1 >= newDeviationNodeIndex);
        assert(deviationNodeIndex <= newDeviationNodeIndex);

        for(size_t i = deviationNodeIndex; i < newDeviationNodeIndex + 1; i++)
        {
            removedNodes[parentPath.path[i]] = true;
        }

        deviationNodeIndex = newDeviationNodeIndex;
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode>(deviationNodeIndex, graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), removedNodes, parentPath);
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getOutEdgesBegin(const NodeType u) const noexcept
    {
//        const auto allEdges = graph.getOutEdges(u);
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), deviationNodeIndex, removedNodes, parentPath);
        else
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), removedNodes);
    }

    template<bool /*fromSourceNode*/, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getOutEdgesEnd(const NodeType u) const noexcept
    {
        return graph.getOutEdgesEnd(u);
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getLightOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode>(deviationNodeIndex, graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), removedNodes, parentPath);
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getLightOutEdgesBegin(const NodeType u) const noexcept
    {
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), deviationNodeIndex, removedNodes, parentPath);
        else
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), removedNodes);
    }

    template<bool /*fromSourceNode*/, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getLightOutEdgesEnd(const NodeType u) const noexcept
    {
        return graph.getLightOutEdgesEnd(u);
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getHeavyOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode>(deviationNodeIndex, graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), removedNodes, parentPath);
    }

    template<bool fromSourceNode, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getHeavyOutEdgesBegin(const NodeType u) const noexcept
    {
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), deviationNodeIndex, removedNodes, parentPath);
        else
            return typename EdgeRange<fromSourceNode>::Iterator(graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), removedNodes);
    }

    template<bool /*fromSourceNode*/, bool /*noExpressEdges*/ = false>
    [[nodiscard]] auto getHeavyOutEdgesEnd(const NodeType u) const noexcept
    {
        return graph.getHeavyOutEdgesEnd(u);
    }

    [[nodiscard]] auto getNumNodes() const noexcept
    {
        return graph.getNumNodes();
    }

    [[nodiscard]] auto getNumEdges() const noexcept
    {
        return graph.getNumEdges();
    }

    [[nodiscard]] auto getPathLength(const std::span<const NodeType> path) const noexcept
    {
        return graph.getPathLength(path);
    }

    [[nodiscard]] WeightType getEdgeWeight(const NodeType u, const NodeType v) const noexcept
    {
        return graph.getEdgeWeight(u, v);
    }

    [[nodiscard]] WeightType getDelta() const noexcept requires DeltaSteppingGraphConcept<GraphType>
    {
        return graph.getDelta();
    }

    [[nodiscard]] WeightType getHeaviestWeight() const noexcept requires DeltaSteppingGraphConcept<GraphType>
    {
        return graph.getHeaviestWeight();
    }
};

#endif //SRC_DYNAMICGRAPH
