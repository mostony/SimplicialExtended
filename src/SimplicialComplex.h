#pragma once

#include "abstract_model.h"
#include <vector>

class SimplicialComplex : public AbstractModel {
public:
  void AddComplex(std::vector<int> complex);

  void RemoveComplex(std::vector<int> complex);

  /// merge two non intersecting SimplicialComplex
  friend void Merge(SimplicialComplex *current, SimplicialComplex *other);

  /// TODO: change method to enum
  static SimplicialComplex *
  CreateCliqueGraph(const std::vector<std::vector<int>> &g, int k,
                    int method = 0, int total_threads = -1);

  std::vector<std::vector<int>> GetMaxSimplices();

  void AddArc(const std::vector<int> &from, const std::vector<int> &to);
};
