#ifndef SRC_TOOLS_STATISTICS
#define SRC_TOOLS_STATISTICS

#include <iostream>
#include <vector>
#include "tools/Output.h"

/**
 * Stores data about a KSP run. It is assumed that it is only used in single thread setting.
 */
template<bool active>
class Statistics;

template<>
class Statistics<false>
{
public:
    explicit Statistics(const unsigned) noexcept
    {}

    template<bool>
    [[nodiscard]] JSON getJSON() const noexcept
    {
        assert(false);
        return {};
    }

    void countSSCShortestDeviation() const noexcept {}

    void countSSCSecondShortestDeviation() const noexcept {}

    void countSSCShortestDeviationLength() const noexcept {}

    void countSSCSecondShortestDeviationLength() const noexcept {}

    void countSsspStoppedEarly() const noexcept {}

    void countTotalPathComputations(const uint64_t) const noexcept {}

    void storeYellowGraphSize(const unsigned, const uint64_t) const noexcept {}

    void storeNumExpressEdges(const unsigned, const uint64_t) noexcept {}

    void storeNumExploredNodes(const unsigned, const uint64_t) const noexcept {}
};

template<>
class Statistics<true>
{
    using CounterType = uint64_t;
    using ListType = std::vector<uint64_t>;

private:
    // SSSP Computations skipped by...
    // ... guessing the shortest deviation
    CounterType SSCShortestDeviation = 0;
    // ... guessing the second shortest deviation
    CounterType SSCSecondShortestDeviation = 0;
    // ... shortest deviation is already to far away
    CounterType SSCShortestDeviationLength = 0;
    // ... second shortest deviation is already to far away
    CounterType SSCSecondShortestDeviationLength = 0;

    // SSSP computation stopped without reaching the target because of kth candidate
    CounterType ssspStoppedEarly = 0;

    CounterType totalPathComputations = 0;

    std::vector<ListType> yellowGraphSizes;
    std::vector<ListType> numExpressEdges;
    std::vector<ListType> exploredNodes;

    static void countGeneric(CounterType& counter, const uint64_t steps = 1) noexcept
    {
        counter += steps;
    }

    static void storeGeneric(ListType& container, const ListType::value_type value) noexcept
    {
        container.push_back(value);
    }

public:
    explicit Statistics([[maybe_unused]] const unsigned k) noexcept
    {
        std::cerr << "WARNING: Statistics collection is on. This could affect performance." << std::endl;

        yellowGraphSizes.resize(k - 1);
        numExpressEdges.resize(k - 1);
        exploredNodes.resize(k - 1);
    }

    template<bool showYellowGraphStats>
    [[nodiscard]] JSON getJSON() const noexcept
    {
        JSON json;
        // SSC = skipped SSSP computation
        json.add("SSCShortestDeviation", SSCShortestDeviation);
        json.add("SSCSecondShortestDeviation", SSCSecondShortestDeviation);
        json.add("SSCShortestDeviationLength", SSCShortestDeviationLength);
        json.add("SSCSecondShortestDeviationLength", SSCSecondShortestDeviationLength);
        // SESC = stopped early SSSP computation
        json.add("ssspStoppedEarly", ssspStoppedEarly);

        json.add("totalPathComputations", totalPathComputations);
        if constexpr(showYellowGraphStats)
        {
            json.add("yellowGraphSizes", yellowGraphSizes);
            json.add("numExpressEdges", numExpressEdges);
        }
        json.add("exploredNodes", exploredNodes);

        return json;
    }

    void countSSCShortestDeviation() noexcept
    {
        countGeneric(SSCShortestDeviation);
    }

    void countSSCSecondShortestDeviation() noexcept
    {
        countGeneric(SSCSecondShortestDeviation);
    }

    void countSSCShortestDeviationLength() noexcept
    {
        countGeneric(SSCShortestDeviationLength);
    }

    void countSSCSecondShortestDeviationLength() noexcept
    {
        countGeneric(SSCSecondShortestDeviationLength);
    }

    void countSsspStoppedEarly() noexcept
    {
        countGeneric(ssspStoppedEarly);
    }

    void countTotalPathComputations(const uint64_t numPaths) noexcept
    {
        countGeneric(totalPathComputations, numPaths);
    }

    void storeYellowGraphSize(const unsigned i, const uint64_t size) noexcept
    {
        storeGeneric(yellowGraphSizes[i], size);
    }

    void storeNumExpressEdges(const unsigned i, const uint64_t size) noexcept
    {
        storeGeneric(numExpressEdges[i], size);
    }

    void storeNumExploredNodes(const unsigned i, const uint64_t numNodes) noexcept
    {
        storeGeneric(exploredNodes[i], numNodes);
    }

    [[nodiscard]] const auto& getYellowGraphSizes() const noexcept
    {
        return yellowGraphSizes;
    }

    [[nodiscard]] const auto& getNumExpressEdges() const noexcept
    {
        return numExpressEdges;
    }

    [[nodiscard]] const auto& getExploredNodes() const noexcept
    {
        return exploredNodes;
    }

    [[nodiscard]] auto getTotalPathComputations() const noexcept
    {
        return totalPathComputations;
    }
};


#endif //SRC_TOOLS_STATISTICS
