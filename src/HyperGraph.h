#pragma once

#include "abstract_model.h"

class HyperGraph : public AbstractModel {
public:
  void AddEdge(const std::vector<int> &edge);

  void RemoveEdge(const std::vector<int> &edge);

  std::vector<std::vector<int>> GetEdges();
};