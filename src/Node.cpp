#include "Node.h"

#include <algorithm>
#include <queue>
#include <stdexcept>
#include <unordered_set>

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

std::vector<Node*> Node::GetAllUpper(int upper_rank) {
    if (upper_rank < rank) {
        throw std::runtime_error("upper rank should be >= current rank");
    }
    std::queue<Node*> queue;
    std::unordered_set<Node*> visited;
    queue.push(this);
    visited.insert(this);
    std::vector<Node*> result;
    while (queue.size()) {
        auto v = queue.front();
        queue.pop();
        if (v->rank == upper_rank) {
            result.push_back(v);
        }
        for (auto u : v->upper) {
            if (u->rank <= upper_rank && !visited.contains(u)) {
                visited.insert(u);
                queue.push(u);
            }
        }
    }
    return std::move(result);
}
