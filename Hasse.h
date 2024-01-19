#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <vector>

#include "Node.h"

class Hasse {
   public:
    void DebugPrintAll() {
        for (auto &[id, node] : mapping_) {
            node.Debug();
        }
    }

    void AddArc(const std::vector<int> &from, const std::vector<int> &to);

    /// remove node and all upper node
    /// that can reach. In other words all upper sets
    void RemoveNode(const std::vector<int> &remove_node);

    void RemoveArc(const std::vector<int> &from, const std::vector<int> &to);

    // add node and all subsets
    void RecursiveAddNode(const std::vector<int> &add_node);

    /// use mapping to get node. If node not in mapping
    /// then create node. Return ptr to Node
    Node &GetNode(const std::vector<int> &node);

   private:
    void RecursiveRemoveNode(const std::vector<int> &node);

    std::map<std::vector<int>, Node> mapping_;
};
