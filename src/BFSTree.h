#ifndef KSP_BSFTREE_H
#define KSP_BSFTREE_H

#include "Graph.h"

#include <queue>


class BFSTree
{

private:
    struct BFSNode
    {
        NodeType depth;
        NodeType predecessor;
        std::vector<NodeType> children;
        NodeType hangingTreeSize;
    };

    std::vector<BFSNode> tree;

public:
    BFSTree() = delete;

    template<GraphConcept GT>
    BFSTree(const GT& graph, const NodeType root)
      : tree(std::vector<BFSNode>(graph.getNumNodes(), {nullNode, nullNode, {}, 0}))
    {
        std::queue<NodeType> queue;
        queue.push(root);
        tree[root].depth = 0;

        while(!queue.empty())
        {
            NodeType currentNode = queue.front();
            queue.pop();

            for(const auto& neighbor : graph.getOutEdges(currentNode))
            {
                const NodeType target = neighbor.getTarget();
                if(tree[target].depth < nullNode)
                    continue;

                set(target, currentNode);
                queue.push(target);
            }
        }

        calculateHangingTreeSizes(root);
    }

    [[nodiscard]] NodeType getDepth(NodeType node) const
    {
        return tree[node].depth;
    }

    [[nodiscard]] NodeType getHangingTreeSize(const NodeType node) const
    {
        return tree[node].hangingTreeSize;
    }

    [[nodiscard]] NodeType getNumChildren(const NodeType node) const
    {
        return static_cast<NodeType>(tree[node].children.size());
    }

private:
    void set(const NodeType node, const NodeType predecessor)
    {
        tree[node].predecessor = predecessor;
        tree[node].depth = tree[predecessor].depth + 1;
        tree[predecessor].children.push_back(node);
    }

    void calculateHangingTreeSizes(const NodeType node)
    {
        NodeType hangingTreeSize = 1;   // the root it self
        for(const auto& child : tree[node].children)
        {
            calculateHangingTreeSizes(child);
            hangingTreeSize += tree[child].hangingTreeSize;
        }

        tree[node].hangingTreeSize = hangingTreeSize;
    }

};


#endif //KSP_BSFTREE_H
