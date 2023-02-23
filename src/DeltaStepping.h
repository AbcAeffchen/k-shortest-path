#ifndef SRC_DELTASTEPPING
#define SRC_DELTASTEPPING

#include "omp.h"
#include <ranges>
#include <valarray>
#include <tuple>
#include <unordered_set>

#include "SsspTree.h"
#include "concepts/graph.h"
#include "concepts/sssp.h"
#include "tools/Output.h"

template<DeltaSteppingGraphConcept GraphType, unsigned numThreads, bool useEarlyStopping, template<typename> typename SSSPTreeType = SSSPTree>
requires SSSPTreeConcept<SSSPTreeType<typename GraphType::WeightType>>
class DeltaStepping;

template<bool useEarlyStopping, DeltaSteppingGraphConcept GraphType, template<typename> typename SSSPTreeType>
using DeltaSteppingSingleThread = DeltaStepping<GraphType, 1, useEarlyStopping, SSSPTreeType>;

template<unsigned numThreads, bool useEarlyStopping, DeltaSteppingGraphConcept GraphType, template<typename> typename SSSPTreeType>
using DeltaSteppingMultiThread = DeltaStepping<GraphType, numThreads, useEarlyStopping, SSSPTreeType>;

template<DeltaSteppingGraphConcept GraphType, unsigned numThreads, bool useEarlyStopping, template<typename> typename SSSPTreeType>
requires SSSPTreeConcept<SSSPTreeType<typename GraphType::WeightType>>
class DeltaStepping
{
    /**
     * The hash function used to distribute nodes to threads.
     * Since we assume a random graph, this simple hash function should work fine.
     * If the graph is not random, a more complex hash function would be needed.
     * Also we assume the number of threads to be a power of two, so the modulus should be fast.
     */
    static constexpr unsigned nodeToThread(const NodeType node) noexcept
    {
        return node % numThreads;
    }

public:
    using DistanceType = typename GraphType::WeightType;
    static constexpr bool earlyStopping = useEarlyStopping;

private:
    enum class COMPUTE_MODE
    {
        FULL, STOP_WHEN_TARGET_IS_FOUND, STOP_WHEN_TARGET_IS_TOO_FAR_AWAY
    };

    /**
     * Used to store the nodes seen in a bucket without duplicates.
     */
    class NodeCache
    {
        // since vector<bool> cannot be used thread safe, we use vector<uint8_t> to store the booleans.
        std::vector<uint8_t> nodeSeen;
        std::array<std::vector<NodeType>, numThreads> nodeSets;

    public:
        NodeCache() = delete;

        explicit NodeCache(const NodeType numNodesInGraph) noexcept
            : nodeSeen(numNodesInGraph, false)
        {
            const auto reserveSize = static_cast<size_t>(static_cast<float>(numNodesInGraph) * 0.05f / static_cast<float>(numThreads)) + 1;

            for(unsigned i = 0; i < numThreads; i++)
            {
                nodeSets[i].reserve(reserveSize);
            }
        }

        void insert(const NodeType node) noexcept
        {
            if(nodeSeen[node])
                return;

            nodeSeen[node] = true;
            nodeSets[static_cast<size_t>(omp_get_thread_num())].push_back(node);
        }

        const auto& getNodes() const noexcept
        {
            return nodeSets[static_cast<size_t>(omp_get_thread_num())];
        }

        [[nodiscard]] bool contains(const NodeType node) const noexcept
        {
            return nodeSeen[node];
        }

        void clear() noexcept
        {
            for(const auto node : nodeSets[static_cast<size_t>(omp_get_thread_num())])
                nodeSeen[node] = false;

            nodeSets[static_cast<size_t>(omp_get_thread_num())].clear();
        }
    };

    /**
     * Used to Buffer distance update requests distributed to threads by destination nodes.
     */
    class RequestBuffer
    {
        std::array<std::array<std::vector<std::tuple<NodeType, NodeType, DistanceType>>, numThreads>, numThreads> buffer;

    public:
        void reserve(size_t size) noexcept
        {
            size /= numThreads * numThreads;
            size++;

            for(size_t i = 0; i < numThreads; i++)
            {
                for(size_t j = 0; j < numThreads; j++)
                {
                    buffer[i][j].reserve(size);
                }
            }
        }

        void insert(const NodeType node, const NodeType predecessor, const DistanceType newDistance) noexcept
        {
            buffer[static_cast<size_t>(omp_get_thread_num())][nodeToThread(node)].emplace_back(node, predecessor, newDistance);
        }

        const auto& getRequests(const unsigned threadId) const noexcept
        {
            return buffer[threadId][static_cast<size_t>(omp_get_thread_num())];
        }

        void clear() noexcept
        {
            const auto threadId = static_cast<size_t>(omp_get_thread_num());
            for(auto& subBuffer : buffer)
            {
                subBuffer[threadId].clear();
            }
        }
    };

    /**
     * The buckets holding nodes having roughly the same tentative distance, depending on delta.
     * The buckets are split into sub buckets for each thread so they can operate independently.
     */
    class BucketList
    {
        class Bucket
        {
            std::array<std::vector<NodeType>, numThreads> threadBuckets;

        public:
            [[nodiscard]] size_t size() const noexcept
            {
                size_t size = 0;
                for(const auto& subBucket : threadBuckets)
                    size += subBucket.size();

                return size;
            }

            [[nodiscard]] bool empty() const noexcept
            {
                for(const auto& subBucket : threadBuckets)
                    if(!subBucket.empty())
                        return false;

                return true;
            }

            void push_back(const NodeType node) noexcept
            {
                threadBuckets[static_cast<size_t>(omp_get_thread_num())].push_back(node);
            }

            const auto& getNodes() const noexcept
            {
                return threadBuckets[static_cast<size_t>(omp_get_thread_num())];
            }

            void clear() noexcept
            {
                threadBuckets[static_cast<size_t>(omp_get_thread_num())].clear();
            }

            void reserve(size_t reserveSize) noexcept
            {
                reserveSize /= numThreads;

                for(auto& bucket : threadBuckets)
                    bucket.reserve(reserveSize);
            }
        };

        const DistanceType deltaReciprocal;
        unsigned currentBucket = 0;       /// the current bucket
        unsigned currentBucketId = 0;     /// the ID of the current bucket, which is `currentBucket % bucketList.size()`.
        std::vector<Bucket> bucketList;

    public:
        BucketList(const NodeType numNodesInGraph, const DistanceType delta, const DistanceType heaviestWeight)
            : deltaReciprocal(static_cast<DistanceType>(1) / delta)
            , bucketList(static_cast<size_t>(heaviestWeight * deltaReciprocal + 1))
        {
            assert(delta > 0);
            assert(heaviestWeight > 0);

            const auto initialBucketSize = static_cast<size_t>(std::pow(numNodesInGraph, 0.25));
            for(auto& bucket : bucketList)
                bucket.reserve(initialBucketSize);
        }

        [[nodiscard]] size_t getCurrentBucketNumber() const noexcept
        {
            return currentBucket;
        }

        /**
         * Advance to the next non-empty bucket.
         * Returns true if such a bucket is found, false otherwise.
         * Make sure that this is only executed by a single thread.
         */
        template<bool advanceAtLeastOneBucket = true>
        void gotoNextNonEmptyBucket() noexcept
        {
            const unsigned currentBucketStart = currentBucketId;

            if constexpr(advanceAtLeastOneBucket)
            {
                currentBucket++; currentBucketId++;
            }

            for(; currentBucketId < bucketList.size(); currentBucket++, currentBucketId++)
            {
                if(!isCurrentBucketEmpty())
                    return;
            }

            for(currentBucketId = 0; currentBucketId < currentBucketStart; currentBucket++, currentBucketId++)
            {
                if(!isCurrentBucketEmpty())
                    return;
            }
        }

        void insert(const NodeType node, const DistanceType tentativeDistance) noexcept
        {
            // find the right bucket
            // depends on delta and the tentative distance
            const auto targetBucket = getBucketNumber(tentativeDistance);
            assert(targetBucket >= currentBucket);
            auto targetBucketId = targetBucket - currentBucket + currentBucketId;
            if(targetBucketId >= bucketList.size())
                targetBucketId -= bucketList.size();

            bucketList[targetBucketId].push_back(node);
        }

        [[nodiscard]] inline size_t getBucketNumber(const DistanceType distance) const noexcept
        {
            return static_cast<size_t>(distance * deltaReciprocal);
        }

        const auto& getNodesFromCurrentBucket() const noexcept
        {
            return bucketList[currentBucketId].getNodes();
        }

        void clearCurrentBucket() noexcept
        {
            bucketList[currentBucketId].clear();
        }

        [[nodiscard]] bool isCurrentBucketEmpty() const noexcept
        {
            return bucketList[currentBucketId].empty();
        }
    };

    const GraphType& graph;

    std::shared_ptr<SSSPTreeType<DistanceType>> ssspTree;

    BucketList bucketList;
    RequestBuffer requestBuffer;
    NodeCache nodeCache;

    NodeType _target = nullNode;
    bool computationFinished = false;

public:
    explicit DeltaStepping(const GraphType& graph) noexcept
        : graph(graph)
        , ssspTree(std::make_shared<SSSPTreeType<DistanceType>>(graph.getNumNodes()))
        , bucketList(graph.getNumNodes(), graph.getDelta(), graph.getHeaviestWeight())
        , nodeCache(graph.getNumNodes())
    {
        requestBuffer.reserve(static_cast<size_t>(std::pow(graph.getNumNodes(), 0.25)));
    }

    /**
     * Computation stops as soon as the target is found or if the target is further away than the maximal allowed distance.
     */
    void compute(const NodeType start, const NodeType target, const DistanceType maxDistance) noexcept
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        _compute<useEarlyStopping ? COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY : COMPUTE_MODE::FULL>(start, target, maxDistance);
    }

    /**
     * Stops as soon as the target is found.
     */
    void compute(const NodeType start, const NodeType target) noexcept
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        _compute<useEarlyStopping ? COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND : COMPUTE_MODE::FULL>(start, target);
    }

    /**
     * Computes the full SSSP tree.
     */
    void compute(const NodeType start) noexcept
    {
        _compute<COMPUTE_MODE::FULL>(start);
        ssspTree->template setChildLists<numThreads>();
    }

    [[nodiscard]] bool pathFound() const noexcept
    {
        assert(_target != nullNode);
        return pathFound(_target);
    }

    [[nodiscard]] bool pathFound(const NodeType target) const noexcept
    {
        assert(_target == nullNode || _target == target);
        return computationFinished && ssspTree->getDistance(target) != Distance<DistanceType>::max;
    }

    [[nodiscard]] DistanceType getDistance() const noexcept
    {
        assert(_target != nullNode);
        return getDistance(_target);
    }

    [[nodiscard]] DistanceType getDistance(const NodeType target) const noexcept
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getDistance(target);
    }

    [[nodiscard]] PathType getPath() const noexcept
    {
        assert(_target != nullNode);
        return getPath(_target);
    }

    [[nodiscard]] PathType getPath(const NodeType target) const noexcept
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getPath(target);
    }

    [[nodiscard]] std::shared_ptr<const SSSPTreeType<DistanceType>> getSSSPTree() const noexcept
    {
        return ssspTree;
    }

private:
    /*
     * How `_compute` works:
     *  - get the current bucket
     *  - for each node in the bucket add update requests into a buffer for each light edge
     *  - relax all update requests in the buffer. This may reinsert nodes in the current bucket.
     *  - repeat until the bucket is empty.
     *  - for each node that was in the bucket add update requests for each heavy edge into a buffer and relax.
     *    This only needs to be done once.
     *  - move to the next bucket and repeat everything until an end condition is met.
     *
     * details:
     *  - the request buffers have numThreads squared number of vectors. so each thread has its own buffer and can assign
     *    requests for the other threads based on the target nodes. So no mutexes are needed.
     *  - whenever nodes get assigned to a request buffer, they are assigned to a thread by a hash function.
     *    The hash function should distribute the nodes close to uniform to the given threads.
     *  - how to remember the nodes that were in the bucket to process the heavy nodes
     *     - don't remember, just directly assign requests for the heavy edges in a separate list. at the end, sort the list and skip if a node was already relaxed.
     *     - store the nodes in a vector, sort and skip over nodes that edges where already added to the request list.
     *     - use a vector and a bitfield of size n to store if a node is already added.
     *     + use an unordered_set (expected constant time), no bitfield to check.
     */
    template<COMPUTE_MODE cm>
    void _compute(const NodeType start, const NodeType target = nullNode, DistanceType maxDistance = Distance<DistanceType>::max)
    {
        size_t finalBucket = maxDistance == Distance<DistanceType>::max
            ? std::numeric_limits<size_t>::max()
            : bucketList.getBucketNumber(maxDistance);

        ssspTree->set(start, 0, nullNode);

        if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
        {
            if(start == target)
            {
                computationFinished = true;
                return;
            }
        }

        for(const auto& e : graph.template getOutEdges<true>(start))
        {
            ssspTree->set(e.getTarget(), e.getWeight(), start);
            bucketList.insert(e.getTarget(), e.getWeight());
        }

        bucketList.template gotoNextNonEmptyBucket<false>();

        // spawn numThreads threads, add barriers to synchronize
#pragma omp parallel num_threads(numThreads) shared(target, maxDistance, finalBucket) default(none)
        {
            /*
             * In order to work as intended, openmp needs to spawn numThread threads.
             * The openmp directive does not guarantee this.
             */
            assert(numThreads == omp_get_num_threads());

            do
            {
                if constexpr(cm == COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY)
                {
                    if(bucketList.getCurrentBucketNumber() > finalBucket)
                        break;
                }

                const auto currentMinDistance = graph.getDelta() * static_cast<DistanceType>(bucketList.getCurrentBucketNumber());

#pragma omp barrier
                // make relax requests for all light edges of nodes in the current bucket
                do
                {
#pragma omp barrier
                    for(const auto node : bucketList.getNodesFromCurrentBucket())
                    {
                        const auto nodeDistance = ssspTree->getDistance(node);

                        // Since this node was added to the bucket its distance was settled via
                        // another edge and can now be ignored
                        if(nodeDistance < currentMinDistance)
                            continue;

                        nodeCache.insert(node);

                        for(const auto& edge: graph.template getLightOutEdges<false>(node))
                        {
                            assert(edge.getWeight() < graph.getDelta());
                            const auto newDistance = nodeDistance + edge.getWeight();
                            // make only requests that can actually improve the distance
                            if(newDistance >= currentMinDistance)
                                requestBuffer.insert(edge.getTarget(), node, newDistance);
                        }
                    }

                    bucketList.clearCurrentBucket();

#pragma omp barrier
                    relaxAllNodesInBuffer();

#pragma omp barrier
                // if there are still nodes in the current bucket, do the same thing again.
                } while(!bucketList.isCurrentBucketEmpty());

                // if the target is settled, we don't need to relax heavy edges.
                if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
                {
                    if(nodeCache.contains(target))
                        break;  // jumps to the end of the parallel block
                }

#pragma omp barrier
                // make relax requests for all heavy edges of nodes in the current bucket
                for(const auto node : nodeCache.getNodes())
                {
                    const auto nodeDistance = ssspTree->getDistance(node);

                    for(const auto& edge : graph.template getHeavyOutEdges<false>(node))
                    {
                        assert(edge.getWeight() >= graph.getDelta());
                        const auto newDistance = nodeDistance + edge.getWeight();
                        // make only requests that are within range
                        if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
                        {
                            if(newDistance >= maxDistance)
                                break;
                        }

                        requestBuffer.insert(edge.getTarget(), node, newDistance);
                    }
                }

                nodeCache.clear();

#pragma omp barrier
                // relax all nodes in the buffer and add all updated nodes into a bucket according to their distance.
                relaxAllNodesInBuffer();

#pragma omp barrier

#pragma omp single
                {
                    if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
                    {
                        const auto targetDist = ssspTree->getDistance(target);
                        if(maxDistance > targetDist)
                        {
                            maxDistance = targetDist;
                            finalBucket = bucketList.getBucketNumber(targetDist);
                        }
                    }

                    bucketList.template gotoNextNonEmptyBucket<true>();
                }

            } while(!bucketList.isCurrentBucketEmpty());
        }   // end of parallel block. All threads will wait here to be joined.

        // during the parallel execution it is not possible to directly return.
        // Thus, a break is used and then the condition is checked again.
        if constexpr(cm == COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY)
        {
            if(bucketList.getCurrentBucketNumber() > finalBucket)
                return;
        }

        computationFinished = true;
    }

    /**
     * Relax all nodes in the buffer and add all updated nodes into a bucket according to their distance.
     * Assumes that each thread calls this function, so each thread can execute all requests
     * it got from the other threads.
     */
    void relaxAllNodesInBuffer() noexcept
    {
        for(unsigned threadId = 0; threadId < numThreads; threadId++)
        {
            for(const auto& [node, predecessor, newDistance] : requestBuffer.getRequests(threadId))
            {
                if(newDistance >= ssspTree->getDistance(node))
                    continue;

                ssspTree->set(node, newDistance, predecessor);
                bucketList.insert(node, newDistance);
            }
        }

        requestBuffer.clear();
    }

};

/**
 * Specialization for a single thread to reduce overhead.
 */
template<DeltaSteppingGraphConcept GraphType, bool useEarlyStopping, template<typename> typename SSSPTreeType>
requires SSSPTreeConcept<SSSPTreeType<typename GraphType::WeightType>>
class DeltaStepping<GraphType, 1, useEarlyStopping, SSSPTreeType>
{
public:
    using DistanceType = typename GraphType::WeightType;
    static constexpr bool earlyStopping = useEarlyStopping;

private:
    enum class COMPUTE_MODE
    {
        FULL, STOP_WHEN_TARGET_IS_FOUND, STOP_WHEN_TARGET_IS_TOO_FAR_AWAY
    };

    /**
     * Used to store the nodes seen in a bucket without duplicates.
     */
    class NodeCache
    {
        std::vector<bool> nodeSeen;
        std::vector<NodeType> nodeSets;

    public:
        NodeCache() = delete;

        explicit NodeCache(const NodeType numNodesInGraph) noexcept
          : nodeSeen(numNodesInGraph, false)
        {
            nodeSets.reserve(static_cast<size_t>(0.2 * numNodesInGraph));
        }

        void insert(const NodeType node) noexcept
        {
            if(nodeSeen[node])
                return;

            nodeSeen[node] = true;
            nodeSets.push_back(node);
        }

        const auto& getNodes() const noexcept
        {
            return nodeSets;
        }

        [[nodiscard]] bool contains(const NodeType node) const noexcept
        {
            return nodeSeen[node];
        }

        void clear() noexcept
        {
            for(const auto node : nodeSets)
                nodeSeen[node] = false;

            nodeSets.clear();
        }
    };

    /**
     * Used to Buffer distance update requests distributed to threads by destination nodes.
     */
    using RequestBuffer = std::vector<std::tuple<NodeType, NodeType, DistanceType>>;

    /**
     * The buckets holding nodes having roughly the same tentative distance, depending on delta.
     * The buckets are split into sub buckets for each thread so they can operate independently.
     */
    class BucketList
    {
        using Bucket = std::vector<NodeType>;

        const DistanceType deltaReciprocal;
        unsigned currentBucket = 0;       /// the current bucket
        unsigned currentBucketId = 0;     /// the ID of the current bucket, which is `currentBucket % bucketList.size()`.
        std::vector<Bucket> bucketList;

    public:
        BucketList(const NodeType numNodesInGraph, const DistanceType delta, const DistanceType heaviestWeight) noexcept
            : deltaReciprocal(static_cast<DistanceType>(1) / delta)
            , bucketList(static_cast<size_t>(heaviestWeight * deltaReciprocal + 1))
        {
            assert(delta > 0);
            assert(heaviestWeight > 0);

            const auto initialBucketSize = static_cast<size_t>(std::pow(numNodesInGraph, 0.25));
            for(auto& bucket : bucketList)
                bucket.reserve(initialBucketSize);
        }

        [[nodiscard]] size_t getCurrentBucketNumber() const noexcept
        {
            return currentBucket;
        }

        /**
         * Advance to the next non-empty bucket.
         * Returns true if such a bucket is found, false otherwise.
         */
        template<bool advanceAtLeastOneBucket = true>
        bool gotoNextNonEmptyBucket() noexcept
        {
            const unsigned currentBucketStart = currentBucketId;

            if constexpr(advanceAtLeastOneBucket)
            {
                currentBucket++; currentBucketId++;
            }

            for(; currentBucketId < bucketList.size(); currentBucket++, currentBucketId++)
            {
                if(!isCurrentBucketEmpty())
                    return true;
            }

            for(currentBucketId = 0; currentBucketId <= currentBucketStart; currentBucket++, currentBucketId++)
            {
                if(!isCurrentBucketEmpty())
                    return true;
            }

            return false;
        }

        void insert(const NodeType node, const DistanceType tentativeDistance) noexcept
        {
            // find the right bucket
            // depends on delta and the tentative distance
            const auto targetBucket = getBucketNumber(tentativeDistance);
            assert(targetBucket >= currentBucket);
            auto targetBucketId = targetBucket - currentBucket + currentBucketId;
            if(targetBucketId >= bucketList.size())
                targetBucketId -= bucketList.size();

            bucketList[targetBucketId].push_back(node);
        }

        [[nodiscard]] inline size_t getBucketNumber(const DistanceType distance) const noexcept
        {
            return static_cast<size_t>(distance * deltaReciprocal);
        }

        const auto& getNodesFromCurrentBucket() const noexcept
        {
            return bucketList[currentBucketId];
        }

        void clearCurrentBucket() noexcept
        {
            bucketList[currentBucketId].clear();
        }

        [[nodiscard]] bool isCurrentBucketEmpty() const noexcept
        {
            return bucketList[currentBucketId].empty();
        }
    };

    const GraphType& graph;

    std::shared_ptr<SSSPTreeType<DistanceType>> ssspTree;

    BucketList bucketList;
    RequestBuffer requestBuffer;
    NodeCache nodeCache;

    NodeType _target = nullNode;
    bool computationFinished = false;

public:
    explicit DeltaStepping(const GraphType& graph) noexcept
        : graph(graph)
        , ssspTree(std::make_shared<SSSPTreeType<DistanceType>>(graph.getNumNodes()))
        , bucketList(graph.getNumNodes(), graph.getDelta(), graph.getHeaviestWeight())
        , nodeCache(graph.getNumNodes())
    {
        requestBuffer.reserve(static_cast<size_t>(std::pow(graph.getNumNodes(), 0.25)));
    }

    /**
     * Computation stops as soon as the target is found or if the target is further away than the maximal allowed distance.
     */
    void compute(const NodeType start, const NodeType target, DistanceType maxDistance) noexcept
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        _compute<useEarlyStopping ? COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY : COMPUTE_MODE::FULL>(start, target, maxDistance);
    }

    /**
     * Stops as soon as the target is found.
     */
    void compute(const NodeType start, const NodeType target) noexcept
    {
        assert(!SSSPTreeType<DistanceType>::twoWayTraversable);
        _target = target;
        _compute<useEarlyStopping ? COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND : COMPUTE_MODE::FULL>(start, target);
    }

    /**
     * Computes the full SSSP tree.
     */
    void compute(const NodeType start) noexcept
    {
        _compute<COMPUTE_MODE::FULL>(start);
        ssspTree->setChildLists();
    }

    [[nodiscard]] bool pathFound() const noexcept
    {
        assert(_target != nullNode);
        return pathFound(_target);
    }

    [[nodiscard]] bool pathFound(const NodeType target) const noexcept
    {
        assert(_target == nullNode || _target == target);
        return computationFinished && ssspTree->getDistance(target) != Distance<DistanceType>::max;
    }

    [[nodiscard]] DistanceType getDistance() const noexcept
    {
        assert(_target != nullNode);
        return getDistance(_target);
    }

    [[nodiscard]] DistanceType getDistance(const NodeType target) const noexcept
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getDistance(target);
    }

    [[nodiscard]] PathType getPath() const noexcept
    {
        assert(_target != nullNode);
        return getPath(_target);
    }

    [[nodiscard]] PathType getPath(const NodeType target) const noexcept
    {
        assert(computationFinished);
        assert(_target == nullNode || _target == target);
        return ssspTree->getPath(target);
    }

    [[nodiscard]] std::shared_ptr<const SSSPTreeType<DistanceType>> getSSSPTree() const noexcept
    {
        return ssspTree;
    }

private:
    /*
     * How `_compute` works:
     *  - get the current bucket
     *  - for each node in the bucket add update requests into a buffer for each light edge
     *  - relax all update requests in the buffer. This may reinsert nodes in the current bucket.
     *  - repeat until the bucket is empty.
     *  - for each node that was in the bucket add update requests for each heavy edge into a buffer and relax.
     *    This only needs to be done once.
     *  - move to the next bucket and repeat everything until an end condition is met.
     */
    template<COMPUTE_MODE cm>
    void _compute(const NodeType start, const NodeType target = nullNode, DistanceType maxDistance = Distance<DistanceType>::max) noexcept
    {
        size_t finalBucket = maxDistance == Distance<DistanceType>::max
            ? std::numeric_limits<size_t>::max()
            : bucketList.getBucketNumber(maxDistance);

        ssspTree->set(start, 0, nullNode);

        if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
        {
            if(start == target)
            {
                computationFinished = true;
                return;
            }
        }

        for(const auto& e : graph.template getOutEdges<true>(start))
        {
            ssspTree->set(e.getTarget(), e.getWeight(), start);
            bucketList.insert(e.getTarget(), e.getWeight());
        }

        bucketList.template gotoNextNonEmptyBucket<false>();

        do
        {
            if constexpr(cm == COMPUTE_MODE::STOP_WHEN_TARGET_IS_TOO_FAR_AWAY)
            {
                if(bucketList.getCurrentBucketNumber() > finalBucket)
                    return;
            }

            const auto currentMinDistance = graph.getDelta() * static_cast<DistanceType>(bucketList.getCurrentBucketNumber());

            do
            {
                // make relax requests for all nodes in the current bucket
                for(const auto node : bucketList.getNodesFromCurrentBucket())
                {
                    const auto nodeDistance = ssspTree->getDistance(node);

                    // Since this node was added to the bucket its distance was settled via
                    // another edge and can now be ignored
                    if(nodeDistance < currentMinDistance)
                        continue;

                    nodeCache.insert(node);

                    for(const auto& edge : graph.template getLightOutEdges<false>(node))
                    {
                        assert(edge.getWeight() < graph.getDelta());
                        const auto newDistance = nodeDistance + edge.getWeight();
                        // make only requests that can actually improve the distance
                        if(newDistance >= currentMinDistance)
                            requestBuffer.emplace_back(edge.getTarget(), node, newDistance);
                    }
                }

                bucketList.clearCurrentBucket();

                relaxAllNodesInBuffer();

            // if there are still nodes in the current bucket, do the same thing again.
            } while(!bucketList.isCurrentBucketEmpty());

            // if the target is settled, we don't need to relax heavy edges.
            if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
            {
                if(nodeCache.contains(target))
                    break;
            }

            // relax heavy edges
            for(const auto node : nodeCache.getNodes())
            {
                const auto nodeDistance = ssspTree->getDistance(node);

                for(const auto& edge : graph.template getHeavyOutEdges<false>(node))
                {
                    assert(edge.getWeight() >= graph.getDelta());
                    const auto newDistance = nodeDistance + edge.getWeight();
                    // make only requests that are within range
                    if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
                    {
                        if(newDistance >= maxDistance)
                            break;
                    }

                    requestBuffer.emplace_back(edge.getTarget(), node, newDistance);
                }
            }

            nodeCache.clear();

            // relax all nodes in the buffer and add all updated nodes into a bucket according to their distance.
            relaxAllNodesInBuffer();

            if constexpr(cm >= COMPUTE_MODE::STOP_WHEN_TARGET_IS_FOUND)
            {
                const auto targetDist = ssspTree->getDistance(target);
                if(maxDistance > targetDist)
                {
                    maxDistance = targetDist;
                    finalBucket = bucketList.getBucketNumber(targetDist);
                }
            }

        } while(bucketList.template gotoNextNonEmptyBucket<true>());

        computationFinished = true;
    }

    /**
     * Relax all nodes in the buffer and add all updated nodes into a bucket according to their distance.
     */
    void relaxAllNodesInBuffer() noexcept
    {
        for(const auto& [node, predecessor, newDistance] : requestBuffer)
        {
            if(newDistance >= ssspTree->getDistance(node))
                continue;

            ssspTree->set(node, newDistance, predecessor);
            bucketList.insert(node, newDistance);
        }

        requestBuffer.clear();
    }

};

#endif //SRC_DELTASTEPPING
