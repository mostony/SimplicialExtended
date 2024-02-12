#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "Node.h"

class Hasse {
   public:
    void DebugPrintAll();

    void AddArc(const std::vector<int> &from, const std::vector<int> &to);

    /// remove node and all upper node
    /// that can reach. In other words all upper sets
    void RemoveNode(const std::vector<int> &remove_node);

    void RemoveArc(const std::vector<int> &from, const std::vector<int> &to);

    // add node and all subsets
    void RecursiveAddNode(const std::vector<int> &add_node);

    /// use mapping to get node. If node not in mapping
    /// then create node. Return ptr to Node
    Node *GetNode(const std::vector<int> &node);

    // TODO: mb UB
    std::vector<Node *> GetMaxFaces();

    int Size();

    friend void Merge(Hasse &current, Hasse &other);

   private:
    void RecursiveRemoveNode(const std::vector<int> &node);
    // mutable std::mutex mtx;
    std::map<std::vector<int>, std::unique_ptr<Node>> mapping_;
};
