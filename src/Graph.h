#pragma once

#include "abstract_model.h"

class Graph : public AbstractModel {
public:
  void AddEdge(int v, int u);

  void RemoveEdge(int v, int u);

  /// TODO: empty graph
  std::vector<std::vector<int>> GetEdges();
};