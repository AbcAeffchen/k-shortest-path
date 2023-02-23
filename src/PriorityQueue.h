#ifndef KSP_PRIORITYQUEUE_H
#define KSP_PRIORITYQUEUE_H

#include <algorithm>
#include <ranges>

#include "concepts/graph.h"


template<WeightConcept WeightType>
class PriorityQueue
{
    using ElementType = std::tuple<WeightType, NodeType, NodeType>;
    std::vector<ElementType> data;

    static constexpr auto cmp = [](const ElementType& lhs, const ElementType& rhs)
        {
            return std::get<0>(lhs) > std::get<0>(rhs);
        };

public:

    PriorityQueue() = default;

    explicit PriorityQueue(const size_t preallocateSize)
    {
        data.reserve(preallocateSize);
    }

    [[nodiscard]] ElementType top() const
    {
        return data.front();
    }

    void pop()
    {
        std::ranges::pop_heap(data, cmp);
        data.pop_back();
    }

    template<bool updateAfterPush = true>
    void push(const ElementType&& element)
    {
        data.push_back(std::move(element));
        if constexpr(updateAfterPush)
            std::ranges::push_heap(data, cmp);
    }

    template<bool updateAfterPush = true>
    void push(const WeightType tentativeDistance, NodeType node, NodeType predecessor)
    {
        data.emplace_back(tentativeDistance, node, predecessor);
        if constexpr(updateAfterPush)
            std::ranges::push_heap(data, cmp);
    }

    void update()
    {
        std::ranges::make_heap(data, cmp);
    }

    [[nodiscard]] bool empty() const
    {
        return data.empty();
    }

    [[nodiscard]] auto size() const
    {
        return data.size();
    }
};



#endif //KSP_PRIORITYQUEUE_H
