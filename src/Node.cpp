#include "Node.h"
#include <algorithm>
#include <stdexcept>

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

void Node::UpdateWeight(double new_weight) {
  if (new_weight <= 0) {
    throw std::runtime_error("weight should be positive");
  }
  weight = new_weight;
}
