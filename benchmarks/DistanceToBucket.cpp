#include <benchmark/benchmark.h>

#include <vector>
#include <cmath>
#include <random>

std::random_device rd;
const auto seed = rd();

template<int version>
size_t test(size_t size, size_t currentBucket, size_t currentBucketId, float tentativeDistance, float deltaReciprocal);

template<>
size_t test<0>(const size_t size, const size_t currentBucket, const size_t currentBucketId, float tentativeDistance, float deltaReciprocal)
{
    const auto targetBucket = static_cast<size_t>(static_cast<double>(tentativeDistance) *
                                                  static_cast<double>(deltaReciprocal));
    auto targetBucketId = (targetBucket - currentBucket + currentBucketId) % size;

    return targetBucketId;
}

template<>
size_t test<1>(const size_t size, const size_t currentBucket, const size_t currentBucketId, float tentativeDistance, float deltaReciprocal)
{
    const auto targetBucket = static_cast<size_t>(std::floor(static_cast<double>(tentativeDistance) *
                                                             static_cast<double>(deltaReciprocal)));
    auto targetBucketId = targetBucket - currentBucket + currentBucketId;
    if(targetBucketId >= size)
        targetBucketId -= size;

    return targetBucketId;
}

template<>
size_t test<2>(const size_t size, const size_t currentBucket, const size_t currentBucketId, float tentativeDistance, float deltaReciprocal)
{
    const auto targetBucket = static_cast<size_t>(static_cast<double>(tentativeDistance) *
                                                  static_cast<double>(deltaReciprocal));
    auto targetBucketId = targetBucket - currentBucket + currentBucketId;
    if(targetBucketId >= size)
        targetBucketId -= size;

    return targetBucketId;
}

template<>
size_t test<3>(const size_t size, const size_t currentBucket, const size_t currentBucketId, float tentativeDistance, float deltaReciprocal)
{
    const auto targetBucket = static_cast<size_t>(tentativeDistance * deltaReciprocal);
    auto targetBucketId = targetBucket - currentBucket + currentBucketId;
    if(targetBucketId >= size)
        targetBucketId -= size;

    return targetBucketId;
}

template<int version>
static void Test(benchmark::State& state)
{
    std::default_random_engine gen(seed);
    std::uniform_int_distribution<size_t> sizeDist(10, 20);
    std::uniform_int_distribution<size_t> bucketDist(0, 1000);
    std::uniform_real_distribution<float> realDist(0, 1);
    std::uniform_real_distribution<float> deltaDist(0.05f, 0.1f);

    auto size = sizeDist(gen);
    auto currentBucket = bucketDist(gen);
    auto currentBucketId = currentBucket % size;
    auto deltaReciprocal = 1.0f / deltaDist(gen);

    for (auto _ : state)
    {
        auto tentativeDistance = realDist(gen);
        benchmark::DoNotOptimize(test<version>(size, currentBucket, currentBucketId, tentativeDistance, deltaReciprocal));
    }
}

BENCHMARK_TEMPLATE(Test, 0);
BENCHMARK_TEMPLATE(Test, 1);
BENCHMARK_TEMPLATE(Test, 2);
BENCHMARK_TEMPLATE(Test, 3);

BENCHMARK_MAIN();