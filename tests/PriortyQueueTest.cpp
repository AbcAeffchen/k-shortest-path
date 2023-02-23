
#include <gtest/gtest.h>

#include "PriorityQueue.h"

#include <random>

class PriorityQueueTest : public ::testing::Test{};

TEST_F(PriorityQueueTest, SinglePush)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<size_t> nodeDist(0, 100);
    std::uniform_real_distribution<double> weightDist(0, 1);
    nodeDist(generator);

    PriorityQueue<double> pq;

    for(int i = 0; i < 100; i++)
        pq.push({weightDist(generator), nodeDist(generator), nodeDist(generator)});

    ASSERT_EQ(pq.size(), 100);

    auto a = pq.top(); pq.pop();
    auto b = pq.top(); pq.pop();

    ASSERT_LE(a, b);

    size_t expectedSize = 98;
    ASSERT_EQ(pq.size(), expectedSize);

    while(!pq.empty())
    {
        a = pq.top();
        pq.pop();
        ASSERT_LE(a, pq.top());
        ASSERT_LE(std::get<0>(a), std::get<0>(pq.top()));
        ASSERT_EQ(pq.size(), --expectedSize);
    }
}

TEST_F(PriorityQueueTest, SequenzPush)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<size_t> nodeDist(0, 100);
    std::uniform_real_distribution<double> weightDist(0, 1);
    nodeDist(generator);

    PriorityQueue<double> pq;

    for(int i = 0; i < 100; i++)
        pq.template push<false>({weightDist(generator), nodeDist(generator), nodeDist(generator)});

    pq.update();

    ASSERT_EQ(pq.size(), 100);

    auto a = pq.top(); pq.pop();
    auto b = pq.top(); pq.pop();

    ASSERT_LE(a, b);

    size_t expectedSize = 98;
    ASSERT_EQ(pq.size(), expectedSize);

    while(!pq.empty())
    {
        a = pq.top();
        pq.pop();
        ASSERT_LE(a, pq.top());
        ASSERT_LE(std::get<0>(a), std::get<0>(pq.top()));
        ASSERT_EQ(pq.size(), --expectedSize);
    }
}