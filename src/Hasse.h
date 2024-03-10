#pragma once

#include <map>
#include <memory>
#include <vector>

#include "Node.h"

class Hasse {
   public:
    void AddArc(const std::vector<int>& from, const std::vector<int>& to);

    /// remove node and all upper nodes
    void RemoveNode(const std::vector<int>& remove_node);

    void RemoveArc(const std::vector<int>& from, const std::vector<int>& to);

    // add node and all lower nodes
    void RecursiveAddNode(const std::vector<int>& add_node);

    /// use mapping to get node. If node not in mapping
    /// then create node. Return ptr to Node
    Node* GetNode(const std::vector<int>& node);

    // TODO: mb UB
    std::vector<Node*> GetMaxFaces();

    /// return number of nodes
    int Size();

    /// merge two non intersecting Hasse diagrams
    friend void Merge(Hasse& current, Hasse& other);

    std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

    std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

   private:
    void RecursiveRemoveNode(const std::vector<int>& node);
    std::map<std::vector<int>, std::unique_ptr<Node>> mapping_;
};
