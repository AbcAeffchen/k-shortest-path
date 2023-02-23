#ifndef KSP_SPECIALGRAPHS_H
#define KSP_SPECIALGRAPHS_H


#include "Graph.h"

template<size_t length>
class Path
{
public:
    using NodeType = uint32_t;
    using EdgeType = Edge<void, true>;
    using EdgeIdType = uint64_t;
    using WeightType = EdgeType::WeightType;
    using DistanceType = WeightType;

    static constexpr bool isWeighted = true;
    static constexpr bool isDirected = true;
    static constexpr bool isDynamic = false;

    [[nodiscard]] NodeType getNumNodes() const noexcept
    {
        return length;
    }

    [[nodiscard]] EdgeIdType getNumEdges() const noexcept
    {
        return length - 1;
    }

    template<bool /*fromSourceNode*/ = false>
    [[nodiscard]] std::vector<EdgeType> getOutEdges(const NodeType node) const noexcept
    {
        if(node >= length)
            return {};

        return { EdgeType(node + 1) };
    }

    [[nodiscard]] DistanceType getPathLength(const PathType& path) const noexcept
    {
        DistanceType len = 0;
        for(unsigned i = 1; i < path.size(); i++)
        {
            assert(path[i] == path[i-1] + 1);
            len++;
        }

        return len;
    }

    [[nodiscard]] WeightType getEdgeWeight(NodeType u, NodeType v) const noexcept
    {
        return 1;
    }

};


#endif //KSP_SPECIALGRAPHS_H
