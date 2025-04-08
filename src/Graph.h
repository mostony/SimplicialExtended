#pragma once

#include "AbstractModel.h"

/// TODO: dimension of graph???
class Graph : public AbstractModel {
   public:
    void AddEdge(int v, int u);

    void RemoveEdge(int v, int u);

    /// TODO: empty graph
    std::vector<std::vector<int>> GetEdges();

    void BuildFromBinary(std::vector<std::vector<int>> binary, bool on_column) override;
};