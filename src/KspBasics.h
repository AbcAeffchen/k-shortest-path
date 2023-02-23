#ifndef KSP_KSPBASICS_H
#define KSP_KSPBASICS_H

#include <vector>
#include <mutex>

#include "Graph.h"


template<typename WeightType>
struct KSPPath
{
    using ForbiddenEdgesListType = std::vector<NodeType>;

    WeightType length;
    unsigned parentPathId;
    NodeType deviationNodeIndex;     // the index of the deviation node in the parent path.
    ForbiddenEdgesListType forbiddenEdges;
    PathType path;

    [[nodiscard]] ForbiddenEdgesListType getExtendedList(NodeType newNode) const noexcept
    {
        ForbiddenEdgesListType extendedForbiddenEdges = forbiddenEdges;
        extendedForbiddenEdges.push_back(newNode);

        return extendedForbiddenEdges;
    }

    auto operator<=>(const KSPPath& rhs) const noexcept
    {
        return std::tie(length, parentPathId, deviationNodeIndex) <=> std::tie(rhs.length, rhs.parentPathId, rhs.deviationNodeIndex) ;
    }
};

template<WeightConcept DistanceType>
class CandidateList
{
    std::vector<KSPPath<DistanceType>> candidates;
    const unsigned maxCandidates;
    unsigned nextCandidate = 0;
    std::mutex pushMutex;

public:
    explicit CandidateList(unsigned maxCandidates) noexcept
      : maxCandidates(maxCandidates)
    {
        assert(maxCandidates > 0);
        candidates.reserve(maxCandidates);
    }

    /**
     * Can be used in multi thread setting.
     */
    void push(const KSPPath<DistanceType>&& candidate) noexcept
    {
        std::lock_guard<std::mutex> lock(pushMutex);

        // check if new candidate is good enough
        if(candidates.size() >= maxCandidates && candidates.back().length < candidate.length)
            return;

        // append candidate at the end or replace the candidate at the back
        if(candidates.size() < maxCandidates)
            candidates.push_back(std::move(candidate));
        else
            candidates[maxCandidates - 1] = std::move(candidate);

        if(candidates.size() < 2) [[unlikely]]
            return;

        // move the new candidate to the right position.
        for(auto i = candidates.size() - 1; i > nextCandidate; i--)
        {
            if(candidates[i - 1] < candidates[i])
                break;

            std::swap(candidates[i - 1], candidates[i]);
        }
    }

    [[nodiscard]] auto pop() noexcept
    {
        return candidates[nextCandidate++];
    }

    [[nodiscard]] bool hasNext() const noexcept
    {
        return nextCandidate < candidates.size();
    }

    [[nodiscard]] DistanceType getPathLengthLimit() const noexcept
    {
        if(candidates.size() >= maxCandidates)
            return candidates.back().length;
        else
            return Distance<DistanceType>::max;
    }
};

template<WeightConcept WeightType>
void operator<<(std::vector<KSPPath<WeightType>>& list, CandidateList<WeightType>& candidates) noexcept
{
    list.push_back(candidates.pop());
}




#endif //KSP_KSPBASICS_H
