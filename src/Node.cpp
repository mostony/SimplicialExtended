#include "Node.h"

Node::Node() {}

Node::Node(const std::vector<int>& node) {
    data = node;
    rank = node.size();
}
