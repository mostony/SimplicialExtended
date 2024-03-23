#include "Node.h"

Node::Node() {}

// TODO: maybe rank = size - 1
// cause simplicial reasons???
Node::Node(const std::vector<int>& node) {
    data = node;
    sort(data.begin(), data.end());
    rank = node.size();
}
