#include "Hasse.h"

#include <queue>
#include <set>

Node *Hasse::GetNode(const std::vector<int> &node) {
    if (!mapping_.count(node)) {
        mapping_[node] = std::unique_ptr<Node>(new Node(node));
    }
    return mapping_[node].get();
}

void Hasse::AddArc(const std::vector<int> &from, const std::vector<int> &to) {
    auto parent = GetNode(from);
    auto son = GetNode(to);
    parent->sons.push_back(to);
    son->parents.push_back(from);
    son->depth = parent->depth + 1;
}

void Hasse::RemoveNode(const std::vector<int> &node) {
    RecursiveRemoveNode(node);
}

void Hasse::RemoveArc(const std::vector<int> &from, const std::vector<int> &to) {
    auto parent = GetNode(from);
    auto son = GetNode(to);
    for (auto it = parent->sons.begin(); it != parent->sons.end(); it++) {
        if (*it == to) {
            parent->sons.erase(it);
            break;
        }
    }

    // remove parent from son list
    for (auto it = son->parents.begin(); it != son->parents.end(); it++) {
        if (*it == from) {
            son->parents.erase(it);
            break;
        }
    }
}

void Hasse::RecursiveRemoveNode(const std::vector<int> &remove_node) {
    std::queue<std::vector<int>> q;
    std::set<std::vector<int>> used;
    q.push(remove_node);
    while (!q.empty()) {
        auto top = q.front();
        q.pop();
        auto node = GetNode(top);

        for (auto nxt : node->sons) {
            if (!used.count(nxt)) {
                used.insert(nxt);
                q.push(nxt);
            }
        }

        for (auto prev : node->parents) {
            RemoveArc(prev, top);
        }

        node->sons.clear();
        node->parents.clear();
    }

    for (auto &node : used) {
        mapping_.erase(node);
    }
}

void Hasse::RecursiveAddNode(const std::vector<int> &add_node) {
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
    for (auto &[id, node] : mapping_) {
        (*node).Debug();
    }
}

std::vector<Node *> Hasse::GetMaxFaces() {
    std::vector<Node *> result;
    for (auto &[id, node] : mapping_) {
        if (node->sons.empty()) {
            result.push_back(node.get());
        }
    }
    return result;
}
