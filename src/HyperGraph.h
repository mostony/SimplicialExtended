#pragma once

#include "abstract_model.h"

class HyperGraph : public AbstractModel {
 public:
  void AddEdge(std::vector<int> edge);

  void RemoveEdge(std::vector<int> edge);

  std::vector<std::vector<int>> GetEdges();
};