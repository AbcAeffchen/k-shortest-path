
#include <gtest/gtest.h>

#include <random>

#include "KspBasics.h"

class CandidateListTest : public ::testing::Test{};

TEST_F(CandidateListTest, LengthSorted)
{
    CandidateList<float> list(5);

    list.push({7.0f, 1, 1, {}, {1,2,3,4}});
    list.push({3.0f, 2, 2, {}, {1,2,3,4}});
    list.push({1.0f, 3, 3, {}, {1,2,3,4}});
    list.push({2.0f, 4, 4, {}, {1,2,3,4}});
    list.push({6.0f, 5, 5, {}, {1,2,3,4}});
    list.push({4.0f, 6, 6, {}, {1,2,3,4}});
    list.push({5.0f, 7, 7, {}, {1,2,3,4}});

    for(int i = 1; i <= 5; i++)
    {
        ASSERT_TRUE(list.hasNext());
        ASSERT_FLOAT_EQ(static_cast<float>(i), list.pop().length);
    }

    ASSERT_FALSE(list.hasNext());
}

TEST_F(CandidateListTest, PathSorted)
{
    CandidateList<float> list(5);

    list.push({1.0f, 7, 1, {}, {1,2,3,4}});
    list.push({1.0f, 3, 2, {}, {1,2,3,4}});
    list.push({1.0f, 1, 3, {}, {1,2,3,4}});
    list.push({1.0f, 2, 4, {}, {1,2,3,4}});
    list.push({1.0f, 6, 5, {}, {1,2,3,4}});
    list.push({1.0f, 4, 6, {}, {1,2,3,4}});
    list.push({1.0f, 5, 7, {}, {1,2,3,4}});

    for(unsigned i = 1; i <= 5; i++)
    {
        ASSERT_TRUE(list.hasNext());
        ASSERT_EQ(i, list.pop().parentPathId);
    }

    ASSERT_FALSE(list.hasNext());
}

TEST_F(CandidateListTest, NodeSorted)
{
    CandidateList<float> list(5);

    list.push({1.0f, 0, 7, {}, {1,2,3,4}});
    list.push({1.0f, 0, 3, {}, {1,2,3,4}});
    list.push({1.0f, 0, 1, {}, {1,2,3,4}});
    list.push({1.0f, 0, 2, {}, {1,2,3,4}});
    list.push({1.0f, 0, 6, {}, {1,2,3,4}});
    list.push({1.0f, 0, 4, {}, {1,2,3,4}});
    list.push({1.0f, 0, 5, {}, {1,2,3,4}});

    for(unsigned i = 1; i <= 5; i++)
    {
        ASSERT_TRUE(list.hasNext());
        ASSERT_EQ(i, list.pop().deviationNodeIndex);
    }

    ASSERT_FALSE(list.hasNext());
}