#include "Node.h"

Node::Node() {}

Node::Node(const std::vector<int>& node) {
    data = node;
    sort(data.begin(), data.end());
    rank = node.size();
}
