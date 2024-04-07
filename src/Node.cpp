#include "Node.h"

Node::Node() {
}

Node::Node(const std::vector<int>& node) {
  data = node;
  sort(data.begin(), data.end());
  rank = static_cast<int>(node.size()) - 1;
  size = node.size();
}

Node::Node(const std::vector<int>& node, int init_rank) {
  data = node;
  sort(data.begin(), data.end());
  rank = init_rank;
  size = node.size();
}
