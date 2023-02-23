
#include <gtest/gtest.h>

#include <chrono>
#include <thread>

#include "tools/Timer.h"

class TimerTest : public ::testing::Test{};

TEST_F(TimerTest, Basic)
{
    const int testTimeMs = 1000;

    Timer timer;
    timer.start();
    std::this_thread::sleep_for(std::chrono::milliseconds(testTimeMs));
    timer.stop();

    auto measuredDifference = std::abs(testTimeMs * 1000 - timer.getTotalDurationUs());

    ASSERT_LE(measuredDifference, 250) << timer << "us";    // 250 microseconds = 0.25 milliseconds
}

TEST_F(TimerTest, Scoped)
{
    const int testTimeMs = 1000;

    Timer timer;

    {
        ScopedTimer sc(timer);
        std::this_thread::sleep_for(std::chrono::milliseconds(testTimeMs));
    }

    auto measuredDifference = std::abs(testTimeMs * 1000 - timer.getTotalDurationUs());

    ASSERT_LE(measuredDifference, 250) << timer << "us";    // 250 microseconds = 0.25 milliseconds
}

TEST_F(TimerTest, MultipleRuns)
{
    const int testTimeMs = 100;
    const int numRuns = 10;
    Timer timer;

    for(int i = 0; i < numRuns; i++)
    {
        ScopedTimer sc(timer);
        std::this_thread::sleep_for(std::chrono::milliseconds(testTimeMs));
    }

    auto measuredDifference = std::abs(testTimeMs * numRuns * 1000 - timer.getTotalDurationUs());

    ASSERT_LE(measuredDifference, 250 * numRuns) << timer << "us";    // 250 microseconds = 0.25 milliseconds
}
