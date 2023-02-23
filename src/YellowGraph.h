#ifndef SRC_YELLOWGRAPH
#define SRC_YELLOWGRAPH

#include <climits>
#include <cmath>

#include "KspBasics.h"
#include "SsspTree.h"
#include "concepts/graph.h"

/**
 * Lightweight wrapper to make a graph dynamic.
 * This wrapper assumes that only very few nodes and edges get removed.
 */
template<GraphConcept GraphType>
class YellowGraph
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

    const GraphType& graph;
    const ParentPathType& parentPath;
    NodeType deviationNodeIndex = 0;
    SSSPTreeTwoWay<DistanceType> const * const reverseSSSPTree;

    // node coloring stuff
    enum class COLOR : uint8_t{green = 1, yellow = 2, red = 3};
    std::vector<COLOR> nodeColor;
    mutable std::unordered_map<NodeType, EdgeType const *> expressEdgeMap;
    mutable std::mutex expressEdgeMapMutex;

public:
    template<bool fromSourceNode, bool noExpressEdges = false>
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
            const YellowGraph<GraphType>& yellowGraph;
            mutable EdgeType expressEdge;
            const NodeType sourceNode;

        public:
            Iterator(const BaseEdgeListIteratorType begin, const BaseEdgeListIteratorType end,
                     const NodeType deviationNodeIndex, const YellowGraph<GraphType>& yellowGraph,
                     const ParentPathType& parentPath, const NodeType sourceNode) noexcept requires fromSourceNode
                : AdditionalIteratorFields(deviationNodeIndex, parentPath)
                , currentEdge(begin), end(end)
                , yellowGraph(yellowGraph)
                , expressEdge(parentPath.path.back(), Weight<WeightType>::max)
                , sourceNode(sourceNode)
            {
                advanceToNextEdgeInRange<false>();
            }

            Iterator(const BaseEdgeListIteratorType& begin, const BaseEdgeListIteratorType& end,
                     const YellowGraph<GraphType>& yellowGraph, const NodeType target,
                     const NodeType sourceNode) noexcept requires (!fromSourceNode)
                : currentEdge(begin), end(end)
                , yellowGraph(yellowGraph)
                , expressEdge(target, 1.0)
                , sourceNode(sourceNode)
            {
                advanceToNextEdgeInRange<false>();
            }

            inline void operator++() noexcept
            {
                advanceToNextEdgeInRange<true>();
            }

        private:
            template<bool atLeastOne>
            void advanceToNextEdgeInRange()
            {
                if constexpr(atLeastOne)
                    ++currentEdge;

                while(currentEdge != end && !currentEdgeIsInRange())
                    ++currentEdge;
            }

        public:
            inline const EdgeType& operator*() const noexcept
            {
                return getEdge();
            }

            inline auto const * operator->() const noexcept
            {
                return &getEdge();
            }

        private:
            const EdgeType& getEdge() const noexcept
            {
                if constexpr(noExpressEdges)
                {
                    return *currentEdge;
                }
                else
                {
                    if(yellowGraph.isGreen(currentEdge->getTarget()))
                    {
                        yellowGraph.updateExpressEdgeMap(sourceNode, &*currentEdge);
                        expressEdge.setWeight(currentEdge->getWeight());
                        return expressEdge;
                    }
                    else
                    {
                        return *currentEdge;
                    }
                }
            }

        public:
            inline friend bool operator!=(const Iterator& iterator, const BaseEdgeListIteratorType& rhs) noexcept
            {
                return iterator.currentEdge != rhs;
            }

            inline friend bool operator==(const Iterator& iterator, const BaseEdgeListIteratorType& rhs) noexcept
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

                if constexpr(noExpressEdges)
                {
                    return !yellowGraph.isRed(currentEdge->getTarget());
                }
                else
                {
                    // skip green nodes that do not improve the express edge and all red nodes
                    return !(yellowGraph.isRed(currentEdge->getTarget())
                             || (yellowGraph.isGreen(currentEdge->getTarget()) && currentEdge->getWeight() >= expressEdge.getWeight()));
                }
           }
        };

        const BaseEdgeListIteratorType allEdgesBegin;
        const BaseEdgeListIteratorType allEdgesEnd;
        const NodeType deviationNodeIndex;
        const YellowGraph<GraphType>& yellowGraph;
        const ParentPathType& parentPath;
        const NodeType sourceNode;

    public:
        EdgeRange(const NodeType deviationNodeIndex,
                  const BaseEdgeListIteratorType allEdgesBegin, const BaseEdgeListIteratorType allEdgesEnd,
                  const YellowGraph<GraphType>& yellowGraph, const ParentPathType& parentPath,
                  const NodeType sourceNode) noexcept
            : allEdgesBegin(allEdgesBegin), allEdgesEnd(allEdgesEnd)
            , deviationNodeIndex(deviationNodeIndex)
            , yellowGraph(yellowGraph)
            , parentPath(parentPath)
            , sourceNode(sourceNode)
        {}

        auto begin() const noexcept
        {
            if constexpr(fromSourceNode)
                return Iterator(allEdgesBegin, allEdgesEnd, deviationNodeIndex, yellowGraph, parentPath, sourceNode);
            else
                return Iterator(allEdgesBegin, allEdgesEnd, yellowGraph, parentPath.path.back(), sourceNode);
        }

        auto end() const noexcept
        {
            return allEdgesEnd;
        }
    };

private:
    [[nodiscard]] COLOR getNodeColor(const NodeType node) const noexcept
    {
        return nodeColor[node];
    }

    void setNodeColor(const NodeType node, const COLOR color) noexcept
    {
        nodeColor[node] = color;
    }

    [[nodiscard]] bool isGreen(const NodeType node) const noexcept
    {
        return nodeColor[node] == COLOR::green;
    }

    [[nodiscard]] bool isRed(const NodeType node) const noexcept
    {
        return nodeColor[node] == COLOR::red;
    }

    /**
     * If SSSP stops as soon as the target is discovered, only a constant number of express
     * edges can be relaxed and need yo be stored.
     */
    void updateExpressEdgeMap(const NodeType node, EdgeType const * const newEdge) const noexcept
    {
        if(!expressEdgeMap.contains(node) || expressEdgeMap.at(node)->getWeight() > newEdge->getWeight())
        {
            std::lock_guard<std::mutex> lock(expressEdgeMapMutex);
            expressEdgeMap[node] = newEdge;
        }
    }

public:
    [[nodiscard]] PathType expressEdgeToPath(const NodeType node) const noexcept
    {
        assert(expressEdgeMap.contains(node));

        return reverseSSSPTree->getReversePath(expressEdgeToOriginalTarget(node));
    }

    [[nodiscard]] NodeType expressEdgeToOriginalTarget(const NodeType node) const noexcept
    {
        assert(expressEdgeMap.contains(node));

        return expressEdgeMap.at(node)->getTarget();
    }

    YellowGraph(const GraphType& graph, const KSPPath<typename GraphType::WeightType>& parentPath,
                SSSPTreeTwoWay<DistanceType> const * const reverseSSSPTree) noexcept
        : graph(graph), parentPath(parentPath), reverseSSSPTree(reverseSSSPTree)
          , nodeColor(graph.getNumNodes(), COLOR::green)
    {
        assert(reverseSSSPTree != nullptr);
    }

    YellowGraph(const GraphType& graph, const KSPPath<typename GraphType::WeightType>& parentPath,
                const unsigned deviationNodeIndex, SSSPTreeTwoWay<DistanceType> const * const reverseSSSPTree) noexcept
        : graph(graph), parentPath(parentPath), reverseSSSPTree(reverseSSSPTree)
        , nodeColor(graph.getNumNodes(), COLOR::green)
    {
        assert(reverseSSSPTree != nullptr);
        advanceDeviationNodeIndex(deviationNodeIndex);
    }

    YellowGraph(const YellowGraph<GraphType>& yellowGraph)
        : graph(yellowGraph.graph), parentPath(yellowGraph.parentPath)
        , deviationNodeIndex(yellowGraph.deviationNodeIndex)
        , reverseSSSPTree(yellowGraph.reverseSSSPTree)
        , nodeColor(yellowGraph.nodeColor)
    {}

    void advanceDeviationNodeIndex(const NodeType newDeviationNodeIndex) noexcept
    {
        assert(parentPath.deviationNodeIndex <= newDeviationNodeIndex && parentPath.path.size() - 1 >= newDeviationNodeIndex);
        assert(deviationNodeIndex <= newDeviationNodeIndex);

        std::vector<NodeType> stack;
        stack.reserve(static_cast<size_t>(std::sqrt(graph.getNumNodes())));
        // color all nodes up to the new deviation node red.
        for(size_t i = deviationNodeIndex; i < newDeviationNodeIndex + 1; i++)
        {
            nodeColor[parentPath.path[i]] = COLOR::red;
        }

        // add all green children of the new red nodes to the stack for coloring.
        for(size_t i = deviationNodeIndex; i < newDeviationNodeIndex + 1; i++)
        {
            for(const auto child : reverseSSSPTree->getChildren(parentPath.path[i]))
            {
                if(nodeColor[child] == COLOR::green)
                {
                    nodeColor[child] = COLOR::yellow;
                    stack.push_back(child);
                }
            }
        }

        deviationNodeIndex = newDeviationNodeIndex;

        findYellowNodes(stack);

        // old express edges are no longer valid
        expressEdgeMap.clear();
    }

    [[nodiscard]] NodeType getNumYellowNodes() const noexcept
    {
        NodeType numYellowNodes = 0;

        for(const auto color : nodeColor)
            if(color == COLOR::yellow)
                numYellowNodes++;

        return numYellowNodes;
    }

private:
    void findYellowNodes(std::vector<NodeType>& stack)
    {
        while(!stack.empty())
        {
            const auto node = stack.back();
            stack.pop_back();

            for(const auto child : reverseSSSPTree->getChildren(node))
            {
                // DFS on the reverse shortest path tree can only encounter green nodes (they become yellow)
                // and red nodes (DFS ignores them).
                if(nodeColor[child] == COLOR::green)
                {
                    nodeColor[child] = COLOR::yellow;
                    stack.push_back(child);
                }
            }
        }
    }

public:
    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode, noExpressEdges>(deviationNodeIndex, graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), *this, parentPath, u);
    }

    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getOutEdgesBegin(const NodeType u) const noexcept
    {
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), deviationNodeIndex, *this, parentPath, u);
        else
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getOutEdgesBegin(u), graph.getOutEdgesEnd(u), *this, u);
    }

    template<bool /*fromSourceNode*/>
    [[nodiscard]] auto getOutEdgesEnd(const NodeType u) const noexcept
    {
        return graph.getOutEdgesEnd(u);
    }

    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getLightOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode, noExpressEdges>(deviationNodeIndex, graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), *this, parentPath, u);
    }

    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getLightOutEdgesBegin(const NodeType u) const noexcept
    {
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), deviationNodeIndex, *this, parentPath, u);
        else
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getLightOutEdgesBegin(u), graph.getLightOutEdgesEnd(u), *this, u);
    }

    template<bool /*fromSourceNode*/>
    [[nodiscard]] auto getLightOutEdgesEnd(const NodeType u) const noexcept
    {
        return graph.getLightOutEdgesEnd(u);
    }

    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getHeavyOutEdges(const NodeType u) const noexcept
    {
        return EdgeRange<fromSourceNode, noExpressEdges>(deviationNodeIndex, graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), *this, parentPath, u);
    }

    template<bool fromSourceNode, bool noExpressEdges = false>
    [[nodiscard]] auto getHeavyOutEdgesBegin(const NodeType u) const noexcept
    {
        if constexpr(fromSourceNode)
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), deviationNodeIndex, *this, parentPath, u);
        else
            return typename EdgeRange<fromSourceNode, noExpressEdges>::Iterator(graph.getHeavyOutEdgesBegin(u), graph.getHeavyOutEdgesEnd(u), *this, u);
    }

    template<bool /*fromSourceNode*/>
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

    [[nodiscard]] size_t getNumExpressEdges() const noexcept
    {
        return expressEdgeMap.size();
    }
};

#endif //SRC_YELLOWGRAPH
