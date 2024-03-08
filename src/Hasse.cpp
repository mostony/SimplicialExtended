#include "Hasse.h"

#include <cassert>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>

Node* Hasse::GetNode(const std::vector<int>& node) {
    if (!mapping_.count(node)) {
        mapping_[node] = std::unique_ptr<Node>(new Node(node));
    }
    return mapping_[node].get();
}

void Hasse::AddArc(const std::vector<int>& from, const std::vector<int>& to) {
    auto parent = GetNode(from);
    auto son = GetNode(to);
    parent->upper.push_back(to);
    son->lower.push_back(from);
    son->depth = parent->depth + 1;
}

void Hasse::RemoveNode(const std::vector<int>& node) {
    RecursiveRemoveNode(node);
}

void Hasse::RemoveArc(const std::vector<int>& from,
                      const std::vector<int>& to) {
    auto parent = GetNode(from);
    auto son = GetNode(to);
    for (auto it = parent->upper.begin(); it != parent->upper.end(); it++) {
        if (*it == to) {
            parent->upper.erase(it);
            break;
        }
    }

    // remove parent from son list
    for (auto it = son->lower.begin(); it != son->lower.end(); it++) {
        if (*it == from) {
            son->lower.erase(it);
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

        for (auto nxt : node->upper) {
            if (!used.count(nxt)) {
                used.insert(nxt);
                q.push(nxt);
            }
        }

        for (auto prev : node->lower) {
            RemoveArc(prev, top);
        }

        node->upper.clear();
        node->lower.clear();
    }

    for (auto& node : used) {
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
            auto lower = GetNode(next);

            AddArc(lower->data, top);

            if (!was_before) {
                q.push(next);
            }
        }
    }
}

void Hasse::DebugPrintAll() {
    for (auto& [id, node] : mapping_) {
        node->Debug();
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
    std::map<std::vector<int>, int> distance;
    q.push(node);
    distance[node] = p;
    std::vector<std::vector<int>> result;
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        int cur_distance = distance[v];
        if (cur_distance == k) {
            result.emplace_back(v);
            continue;
        }
        for (auto nxt : GetNode(v)->upper) {
            if (!distance.count(nxt)) {
                distance[nxt] = cur_distance + 1;
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
    std::map<std::vector<int>, int> distance;
    for (auto source : sources) {
        assert(source.size() == k);
        q.push(source);
        distance[source] = k;
    }
    std::vector<std::vector<int>> result;
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        int cur_distance = distance[v];
        if (cur_distance == p) {
            result.emplace_back(v);
            continue;
        }
        for (auto nxt : GetNode(v)->lower) {
            if (!distance.count(nxt)) {
                distance[nxt] = cur_distance - 1;
                q.push(nxt);
            }
        }
    }
    return result;
}