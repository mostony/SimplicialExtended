#pragma once

#include <vector>

#include "abstract_model.h"

class CombinatorialComplex : public AbstractModel {
 public:
  void Build(std::vector<std::vector<int>> data);

  std::vector<std::vector<int>> GetSubsets();
};
