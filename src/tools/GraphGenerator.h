#ifndef KSP_GENERATOR_H
#define KSP_GENERATOR_H

#include "Graph.h"

#include <algorithm>
#include <iterator>
#include <random>
#include <stack>
#include <tuple>
#include "tools/Output.h"

class GraphGenerator
{
public:
    using EdgeType = Edge<float, true>;

    using AdjacencyListType = std::vector<std::vector<EdgeType>>;

    static void printAdjacencyList(const AdjacencyListType& al) noexcept
    {
        for(NodeType i = 0; i < al.size(); i++)
        {
            std::cout << "\n" << i << " -> ";
            for(const auto& target : al[i])
                std::cout << target.getTarget() << " (" << target.getWeight() << "), ";
        }
        std::cout << "\n\n";
    }

    static AdjacencyListType gilbert(const bool weighted, const bool directed, NodeType numNodes, double avgDegree) noexcept
    {
        std::random_device seedGen;
        unsigned long seed = seedGen();
        return gilbert(weighted, directed, numNodes, avgDegree, seed);
    }

    /**
     * Generates a graph G(V,E) without self- or multi-edges following the G(n,p) model with p = avgDegree / numNodes.
     * Optimized for sparse graphs running in O(|E|).
     *
     * @tparam weighted if true, weights are drawn uniformly at random from [0,1).
     * @tparam directed
     * @param numNodes
     * @param avgDegree out-degree on directed graphs.
     * @param seed If zero the random_device is used to get a seed.
     * @return
     */
    static AdjacencyListType gilbert(const bool weighted, const double directness, NodeType numNodes, double avgDegree, const unsigned long seed) noexcept
    {
        assert(directness >= 0.0 && directness <= 1.0);

        const bool strictDirected = directness >= 1.0;
        const bool strictUndirected = directness <= 0.0;

        AdjacencyListType neighbors(numNodes);

        const double p_target = avgDegree / numNodes;
        const double p = p_target / (2.0 - directness);

        std::default_random_engine gen(seed);
        std::bernoulli_distribution addReverseEdge(1.0 - directness);
        std::geometric_distribution<EdgeIdType> nextEdge(p);
        std::uniform_real_distribution<float> newWeight(0, 1);

        NodeType u = 0;
        EdgeIdType v = nextEdge(gen);

        normalizeSkipNodes(u, v, numNodes);

        size_t numUndirectedEdges = 0;
        while(u < numNodes)
        {
            float weight = weighted ? newWeight(gen) : 1.0f;

            neighbors[u].emplace_back(static_cast<NodeType>(v), weight);

            if(strictUndirected || (!strictDirected && addReverseEdge(gen)))
            {
                neighbors[v].emplace_back(u, weight);
                numUndirectedEdges++;
            }
            v += nextEdge(gen) + 1;
            normalizeSkipNodes(u, v, numNodes);
        }

        if(strictDirected)
            removeMultiEdges<false>(neighbors);
        else
            removeMultiEdges<true>(neighbors);

        Print::info() << "Undirected: " << numUndirectedEdges << " (" << (numUndirectedEdges * 2) << ")" << std::endl;

//        if(strictUndirected)
//            getGiantComponentUndirected(neighbors);   // todo implement the faster version for undirected graphs?
//        else
        getGiantComponentDirected(neighbors);

        return neighbors;
    }

    static AdjacencyListType grid(const bool weighted, const bool directed, const NodeType width, const NodeType height, const double edgeProb) noexcept
    {
        std::random_device seedGen;
        unsigned long seed = seedGen();
        return grid(weighted, directed, width, height, edgeProb, seed);
    }

    static AdjacencyListType grid(const bool weighted, const bool directed, const NodeType width, const NodeType height, const double edgeProb, const unsigned long seed) noexcept
    {
        std::default_random_engine gen(seed);
        std::bernoulli_distribution addEdge(edgeProb);
        std::uniform_real_distribution<float> newWeight(0, 1);

        const size_t numNodes = static_cast<size_t>(width) * static_cast<size_t>(height);
        assert(numNodes < (1ul << 32));

        // create node map to get random node IDs while still having a grid
        std::vector<NodeType> nodeMap(numNodes, 0);
        for(NodeType i = 1; i < numNodes; i++)
            nodeMap[i] = i;

        std::ranges::shuffle(nodeMap, gen);

        auto coordToId = [width](const NodeType w, const NodeType h){ return h * width + w; };

        AdjacencyListType neighbors(numNodes);
        for(NodeType h = 0; h < height; h++)
        {
            for(NodeType w = 0; w < width; w++)
            {
                NodeType currentNode = nodeMap[coordToId(w, h)];

                // edge to the west
                if(directed && w != 0 && addEdge(gen))
                    storeEdge(neighbors, currentNode, nodeMap[coordToId(w - 1, h)],
                              weighted ? newWeight(gen) : 1.0f, directed);

                // edge to the east
                if(w != width - 1 && addEdge(gen))
                    storeEdge(neighbors, currentNode, nodeMap[coordToId(w + 1, h)],
                              weighted ? newWeight(gen) : 1.0f, directed);

                // edge to the south
                if(directed && h != 0 && addEdge(gen))
                    storeEdge(neighbors, currentNode, nodeMap[coordToId(w, h - 1)],
                              weighted ? newWeight(gen) : 1.0f, directed);

                // edge to the north
                if(h != height - 1 && addEdge(gen))
                    storeEdge(neighbors, currentNode, nodeMap[coordToId(w, h + 1)],
                              weighted ? newWeight(gen) : 1.0f, directed);
            }
        }

        getGiantComponentDirected(neighbors);

        return neighbors;
    }

private:

    static void storeEdge(AdjacencyListType& neighbors, const NodeType source, const NodeType target, const float weight, const bool directed) noexcept
    {
        neighbors[source].emplace_back(target, weight);
        if(!directed)
        {
            neighbors[target].emplace_back(source, weight);
        }
    }

private:
    static void normalizeSkipNodes(NodeType& u, EdgeIdType& v, const NodeType numNodes) noexcept
    {
        while(v >= numNodes)  // compare SkipNeighborBenchmark
        {
            u++;
            v -= numNodes;
        }

        if(u == v)
        {
            if(v == numNodes - 1)
            {
                u++;
                v = 0;
            }
            else
            {
                v++;
            }
        }
    }

    /**
     * Removes multi edges from the adjacency list. Checks if a multi edge is a undirected edge and if so removes both.
     * Returns the number of edges removed.
     */
    template<bool checkForReverseEdges>
    static void removeMultiEdges(AdjacencyListType& adjacencyList) noexcept
    {
        EdgeIdType numRemovedEdges = 0;

        for(size_t u = 0; u < adjacencyList.size(); u++)
        {
            std::ranges::sort(adjacencyList[u], [](const EdgeType& l, const EdgeType& r)
            {
                return l.getTarget() < r.getTarget();
            });

            size_t newRemovedNodes = 0;

            if constexpr(checkForReverseEdges)
            {
                for(size_t i = 1; i < adjacencyList[u].size(); i++)
                {
                    if(adjacencyList[u][i-1].getTarget() != adjacencyList[u][i].getTarget())
                        continue;

                    auto& edge = adjacencyList[u][i];
                    newRemovedNodes += std::erase_if(adjacencyList[edge.getTarget()], [u, edge](const EdgeType& e){
                        return e.getTarget() == u && abs(e.getWeight() - edge.getWeight()) < std::numeric_limits<float>::epsilon();
                    });
                }
            }

            const auto [first, last] = std::ranges::unique(adjacencyList[u], [](const EdgeType& l, const EdgeType& r)
            {
                return l.getTarget() == r.getTarget();
            });

            newRemovedNodes += static_cast<size_t>(std::distance(first, last));
            numRemovedEdges += newRemovedNodes;

            adjacencyList[u].erase(first, last);
        }

        Print::info() << "Multi edges removed: " << numRemovedEdges << std::endl;
    }

    struct NodeProperties{
        NodeType currentComponentId = nullNode;
        NodeType lowLink = nullNode;
        bool onStack = false;
        bool visited = false;
    };

    using NodePropertyList = std::vector<NodeProperties>;

public:
    /**
     * Uses Tarjan's algorithm to find the largest strongly connected component.
     */
    static void getGiantComponentDirected(AdjacencyListType& adjacencyList) noexcept
    {
        std::stack<NodeType> nodeStack;
        NodeType currentComponentId = 0;
        NodePropertyList nodeProperties(adjacencyList.size());
        // (size, id) pairs
        std::deque<std::pair<NodeType, NodeType>> componentSizes;

        for(NodeType v = 0; v < adjacencyList.size(); v++)
        {
            if(nodeProperties[v].currentComponentId == nullNode)
                tarjanIterative(componentSizes, nodeStack, nodeProperties, adjacencyList, v, currentComponentId);
        }

        const auto maxComponentId = std::max_element(componentSizes.begin(), componentSizes.end())->second;

        std::vector<NodeType> nodesToRemove;
        nodesToRemove.reserve(adjacencyList.size() / 20);   // about 5%
        for(NodeType u = 0; u < adjacencyList.size(); u++)
        {
            if(nodeProperties[u].lowLink != maxComponentId)
                nodesToRemove.push_back(u);
        }

        filterGraph(adjacencyList, nodesToRemove);
    }

private:
    static void tarjanIterative(std::deque<std::pair<NodeType, NodeType>>& componentSizes,
                                std::stack<NodeType>& nodeStack,
                                NodePropertyList& nodeProperties,
                                const AdjacencyListType& adjacencyList,
                                const NodeType start, NodeType& currentComponentId) noexcept
    {
        std::stack<std::tuple<NodeType, NodeType>> callStack;
        callStack.push({start, start});

        while(!callStack.empty())
        {
            auto [node, predecessor] = callStack.top();
            nodeProperties[node].visited = true;

            if(!nodeProperties[node].onStack)
            {
                nodeStack.push(node);
                nodeProperties[node].onStack = true;
            }

            if(nodeProperties[node].currentComponentId == nullNode)
            {
                nodeProperties[node].currentComponentId = currentComponentId;
                nodeProperties[node].lowLink = currentComponentId;
                currentComponentId++;
            }

            for(const auto& neighbor: adjacencyList[node])
            {
                if(neighbor.getTarget() == predecessor)
                    continue;

                if(!nodeProperties[neighbor.getTarget()].visited)
                {
                    callStack.push({neighbor.getTarget(), node});
                    break;
                }

                if(nodeProperties[neighbor.getTarget()].onStack)
                    nodeProperties[node].lowLink = std::min(nodeProperties[node].lowLink,
                                                            nodeProperties[neighbor.getTarget()].lowLink);
            }

            if(std::get<0>(callStack.top()) != node)
                continue;

            if(nodeProperties[node].lowLink == nodeProperties[node].currentComponentId)
            {
                NodeType componentSize = 0;
                while(nodeProperties[node].onStack)
                {
                    nodeProperties[nodeStack.top()].lowLink = nodeProperties[node].lowLink;
                    nodeProperties[nodeStack.top()].onStack = false;
                    nodeStack.pop();
                    componentSize++;
                }

                componentSizes.emplace_back(componentSize, nodeProperties[node].lowLink);
            }

            callStack.pop();
        }
    }

public:
    /**
     * Removes nodes from the graph and renumbers the remaining nodes such that they are consecutive again.
     * @param adjacencyList
     * @param nodesToRemove Sorted list of nodes. It is assumed that this list is short compared
     * to the total number of nodes.
     */
    static void filterGraph(AdjacencyListType& adjacencyList, std::vector<NodeType>& nodesToRemove) noexcept
    {
        if(nodesToRemove.empty())
            return;

        std::vector<NodeType> newNodeIds(adjacencyList.size(), nullNode);
        // todo sorting can be skipped if it is ensured that the list already comes sorted.
        std::sort(nodesToRemove.begin(), nodesToRemove.end());

        // move neighbors of removed nodes to the end.
        NodeType currentNewNodeId = 0;
        for( ; currentNewNodeId < nodesToRemove[0]; currentNewNodeId++)
            newNodeIds[currentNewNodeId] = currentNewNodeId;

        for(size_t i = 1; i < nodesToRemove.size(); i++)
        {
            for(size_t j = nodesToRemove[i-1] + 1; j < nodesToRemove[i]; j++, currentNewNodeId++)
            {

                std::swap(adjacencyList[currentNewNodeId], adjacencyList[j]);
                newNodeIds[j] = currentNewNodeId;
            }
        }

        for(size_t j = nodesToRemove[nodesToRemove.size()-1] + 1; j < adjacencyList.size(); j++)
        {
            std::swap(adjacencyList[currentNewNodeId], adjacencyList[j]);
            newNodeIds[j] = currentNewNodeId++;
        }

        // remove all nodes from the list that were moved to the end.
        adjacencyList.erase(adjacencyList.begin() + currentNewNodeId, adjacencyList.end());

        // fix neighbor IDs and remove neighbors that are no longer around.
        for(auto& neighbors : adjacencyList)
        {
            // map old to new IDs
            for(auto& u : neighbors)
            {
                u.setTarget(newNodeIds[u.getTarget()]);
            }

            // remove neighbors that no longer exist
            for(auto it = neighbors.begin(); it != neighbors.end(); )
            {
                if(it->getTarget() == nullNode)
                    neighbors.erase(it);    // num neighbors are assumed to be small.
                else
                    it++;
            }
        }
    }

};


#endif //KSP_GENERATOR_H
