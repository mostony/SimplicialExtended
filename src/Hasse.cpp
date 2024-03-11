#include "Hasse.h"

#include <cassert>
#include <queue>
#include <set>
#include <stdexcept>

Node* Hasse::GetNode(const std::vector<int>& node) {
    if (!mapping_.count(node)) {
        mapping_[node] = std::unique_ptr<Node>(new Node(node));
    }
    return mapping_[node].get();
}

void Hasse::AddArc(const std::vector<int>& from, const std::vector<int>& to) {
    auto low = GetNode(from);
    auto up = GetNode(to);
    low->upper.push_back(to);
    up->lower.push_back(from);
}

void Hasse::RemoveNode(const std::vector<int>& node) {
    RecursiveRemoveNode(node);
}

void Hasse::RemoveArc(const std::vector<int>& from,
                      const std::vector<int>& to) {
    auto low = GetNode(from);
    auto up = GetNode(to);
    for (auto it = low->upper.begin(); it != low->upper.end(); it++) {
        if (*it == to) {
            low->upper.erase(it);
            break;
        }
    }

    for (auto it = up->lower.begin(); it != up->lower.end(); it++) {
        if (*it == from) {
            up->lower.erase(it);
            break;
        }
    }
}

void Hasse::RecursiveRemoveNode(const std::vector<int>& remove_node) {
    std::queue<std::vector<int>> q;
    std::set<std::vector<int>> used;
    q.push(remove_node);
    while (!q.empty()) {
        auto top = q.front();
        q.pop();
        auto node = GetNode(top);

        for (const auto& nxt : node->upper) {
            if (!used.count(nxt)) {
                used.insert(nxt);
                q.push(nxt);
            }
        }

        for (const auto& prev : node->lower) {
            RemoveArc(prev, top);
        }

        node->upper.clear();
        node->lower.clear();
    }

    for (const auto& node : used) {
        mapping_.erase(node);
    }
}

void Hasse::RecursiveAddNode(const std::vector<int>& add_node) {
    if (mapping_.count(add_node)) {
        return;
    }

    std::queue<std::vector<int>> q;
    q.push(add_node);
    while (!q.empty()) {
        auto top = q.front();
        q.pop();

        // node is leaf
        if (top.size() == 1) {
            continue;
        }
        auto node = GetNode(top);

        for (size_t i = 0; i < top.size(); i++) {
            auto next = top;
            next.erase(next.begin() + i);

            bool was_before = mapping_.count(next);
            auto low = GetNode(next);

            AddArc(low->data, top);

            if (!was_before) {
                q.push(next);
            }
        }
    }
}

std::vector<Node*> Hasse::GetMaxFaces() {
    std::vector<Node*> result;
    for (auto& [id, node] : mapping_) {
        if (node->upper.empty()) {
            result.push_back(node.get());
        }
    }
    return result;
}

int Hasse::Size() {
    return mapping_.size();
}

void Merge(Hasse& current, Hasse& other) {
    for (auto& [key, value] : other.mapping_) {
        current.mapping_[key] = std::move(value);
    }
}

std::vector<std::vector<int>> Hasse::Incidence(std::vector<int> node, int k) {
    int p = node.size();
    if (k < p) {
        throw std::runtime_error("k < size of node");
    }
    std::queue<std::vector<int>> q;
    std::set<std::vector<int>> visited;

    q.push(node);
    visited.insert(node);
    std::vector<std::vector<int>> result;
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        if (GetNode(v)->rank == k) {
            result.emplace_back(v);
            continue;
        }
        for (auto nxt : GetNode(v)->upper) {
            if (!visited.count(nxt)) {
                visited.insert(nxt);
                q.push(nxt);
            }
        }
    }
    return result;
}

std::vector<std::vector<int>> Hasse::Degree(std::vector<int> node, int k) {
    int p = node.size();
    auto sources = Incidence(node, k);
    std::queue<std::vector<int>> q;
    std::set<std::vector<int>> visited;
    for (auto source : sources) {
        assert(source.size() == k);
        q.push(source);
        visited.insert(source);
    }
    std::vector<std::vector<int>> result;
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        if (GetNode(v)->rank == p) {
            result.emplace_back(v);
            continue;
        }
        for (auto nxt : GetNode(v)->lower) {
            if (!visited.count(nxt)) {
                visited.insert(nxt);
                q.push(nxt);
            }
        }
    }
    return result;
}