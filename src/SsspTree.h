#ifndef KSP_SSSPTREE_H
#define KSP_SSSPTREE_H

#include "Graph.h"
#include "omp.h"
#include <algorithm>

template<WeightConcept DistanceType>
struct SSSPNodeWithAnalytics
{
    DistanceType distance = Distance<DistanceType>::max;
    NodeType predecessor = nullNode;
    NodeType depth = nullNode;
    NodeType hangingTreeSize = 0;
    std::vector<NodeType> children = {};
};

template<WeightConcept DistanceType>
struct SSSPNodeTwoWay
{
    DistanceType distance = Distance<DistanceType>::max;
    NodeType predecessor = nullNode;
    std::vector<NodeType> children = {};
};

template<WeightConcept DistanceType>
struct SSSPNodeKIM
{
    DistanceType distance = Distance<DistanceType>::max;
    NodeType predecessor = nullNode;
    NodeType deviationIndex = nullNode;
    std::vector<NodeType> children = {};
};

template<WeightConcept DistanceType>
struct SSSPNodeBasic
{
    DistanceType distance = Distance<DistanceType>::max;
    NodeType predecessor = nullNode;
};

template<WeightConcept DT, template<WeightConcept> typename SSSPNodeTemplate>
class SSSPTreeTemplate
{
public:
    using DistanceType = DT;

private:
    using SSSPNodeType = SSSPNodeTemplate<DistanceType>;

public:
    static constexpr bool twoWayTraversable = std::is_same_v<SSSPNodeType, SSSPNodeTwoWay<DistanceType>> ||
                                              std::is_same_v<SSSPNodeType, SSSPNodeKIM<DistanceType>> ||
                                              std::is_same_v<SSSPNodeType, SSSPNodeWithAnalytics<DistanceType>>;
    static constexpr bool containsAnalytics = std::is_same_v<SSSPNodeType, SSSPNodeWithAnalytics<DistanceType>>;
    static constexpr bool kimReady          = std::is_same_v<SSSPNodeType, SSSPNodeKIM<DistanceType>>;

private:
    std::vector<SSSPNodeType> tree;

public:
    SSSPTreeTemplate() = delete;

    explicit SSSPTreeTemplate(NodeType size) noexcept
      : tree(std::vector<SSSPNodeType>(size))
    {}

    void set(const NodeType node, const DistanceType distance) noexcept
    {
        tree[node].distance = distance;
    }

    void set(const NodeType node, const DistanceType distance, const NodeType predecessor) noexcept
    {
        set(node, distance);
        tree[node].predecessor = predecessor;
    }

    void computeNodeDepth()
    {
        // todo
        assert(false);
    }

    template<unsigned numThreads = 1>
    void setChildLists() noexcept
    {
        if constexpr (twoWayTraversable)
        {
            if constexpr(numThreads == 1)
            {
                std::vector<std::tuple<NodeType, NodeType>> buffer;
                buffer.reserve(tree.size());

                for(NodeType i = 0; i < tree.size(); i++)
                {
                    if(tree[i].predecessor != nullNode)
                        buffer.emplace_back(i, tree[i].predecessor);
                }

                std::ranges::sort(buffer, [](const auto& l, const auto& r) {
                    return std::get<1>(l) < std::get<1>(r);
                });

                for(const auto&[node, predecessor]: buffer)
                {
                    tree[predecessor].children.push_back(node);
                }
            }
            else
            {
                const double chunkRatio = static_cast<double>(numThreads) / static_cast<double>(tree.size());
                std::array<std::vector<std::tuple<NodeType, NodeType>>, numThreads * numThreads> buffer;
                for(auto & b : buffer)
                    b.reserve(tree.size() / (numThreads * numThreads));

#pragma omp parallel num_threads(numThreads) shared(chunkRatio, buffer) default(none)
                {
#pragma omp for
                    for(NodeType i = 0; i < tree.size(); i++)
                    {
                        if(tree[i].predecessor != nullNode)
                            buffer[getBufferId<numThreads>(omp_get_thread_num(), static_cast<double>(tree[i].predecessor) * chunkRatio)]
                                .emplace_back(i, tree[i].predecessor);
                    }

                    for(unsigned chunkId = 0; chunkId < numThreads; chunkId++)
                    {
                        std::sort(buffer[getBufferId<numThreads>(chunkId, omp_get_thread_num())].begin(),
                                  buffer[getBufferId<numThreads>(chunkId, omp_get_thread_num())].end(),
                                  [](const auto& l, const auto& r) {
                                    return std::get<1>(l) < std::get<1>(r);
                        });

                        for(const auto&[node, predecessor]: buffer[getBufferId<numThreads>(chunkId, omp_get_thread_num())])
                        {
                            tree[predecessor].children.push_back(node);
                        }
                    }
                }
            }
        }
    }

private:
    template<unsigned numThreads, typename T, typename S>
    requires (std::is_integral_v<T> || std::is_floating_point_v<T>)
        && (std::is_integral_v<S> || std::is_floating_point_v<S>)
    static constexpr unsigned getBufferId(T x, S y) noexcept
    {
        return static_cast<unsigned>(x) * numThreads + static_cast<unsigned>(y);
    }

public:

    [[nodiscard]] DistanceType getDistance(const NodeType node) const noexcept
    {
        return tree[node].distance;
    }

    [[nodiscard]] PathType getReversePath(NodeType node) const noexcept
    {
        if(tree[node].predecessor == nullNode)
            return {node};

        PathType reversePath;     // todo reserve some space? maybe something like log n elements?

        do
        {
            reversePath.push_back(node);
            node = tree[node].predecessor;
        }
        while(node != nullNode);

        // the tree stores for each node its predecessor. So the path is reversed by default.

        return reversePath;
    }

    [[nodiscard]] PathType getPath(NodeType node) const noexcept
    {
        auto reversePath = getReversePath(node);

        std::ranges::reverse(reversePath);

        return reversePath;
    }

    void calculateHangingTreeSizes(NodeType node) noexcept requires requires(SSSPNodeType x){ x.hangingTreeSize; x.children; }
    {
        NodeType hangingTreeSize = 1;   // the root itself
        for(const auto& child : tree[node].children)
        {
            calculateHangingTreeSizes(child);
            hangingTreeSize += tree[child].hangingTreeSize;
        }

        tree[node].hangingTreeSize = hangingTreeSize;
    }

    [[nodiscard]] NodeType getDepth(const NodeType node) const noexcept requires containsAnalytics
    {
        return tree[node].depth;
    }

    [[nodiscard]] NodeType getHangingTreeSize(const NodeType node) const noexcept requires containsAnalytics
    {
        return tree[node].hangingTreeSize;
    }

    [[nodiscard]] NodeType getNumChildren(const NodeType node) const noexcept requires twoWayTraversable
    {
        return static_cast<NodeType>(tree[node].children.size());
    }

    [[nodiscard]] const auto& getChildren(const NodeType node) const noexcept requires twoWayTraversable
    {
        return tree[node].children;
    }

    [[nodiscard]] auto getPredecessor(const NodeType node) const noexcept
    {
        return tree[node].predecessor;
    }

    [[nodiscard]] auto getDeviationIndex(const NodeType node) const noexcept requires kimReady
    {
        return tree[node].deviationIndex;
    }

    [[nodiscard]] NodeType getNumDiscoveredNodes() const noexcept
    {
        NodeType numDiscoveredNodes = 0;

        for(const auto& node : tree)
            if(node.distance != Distance<DistanceType>::max)
                numDiscoveredNodes++;

        return numDiscoveredNodes;
    }
};

template<WeightConcept DistanceType>
using SSSPTree = SSSPTreeTemplate<DistanceType, SSSPNodeBasic>;

template<WeightConcept DistanceType>
using SSSPTreeTwoWay = SSSPTreeTemplate<DistanceType, SSSPNodeTwoWay>;

template<WeightConcept DistanceType>
using SSSPTreeWithAnalytics = SSSPTreeTemplate<DistanceType, SSSPNodeWithAnalytics>;


#endif //KSP_SSSPTREE_H
