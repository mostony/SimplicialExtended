#pragma once

#include "abstract_model.h"
#include <vector>

class SimplicialComplex : public AbstractModel {
 public:
  void AddSimplex(std::vector<int> simplex);

  void RemoveSimplex(std::vector<int> simplex);

  /// merge two non intersecting SimplicialComplex
  friend void Merge(SimplicialComplex* current, SimplicialComplex* other);

  /// TODO: change method to enum
  static SimplicialComplex* CreateCliqueGraph(
      const std::vector<std::vector<int>>& g, int k, int method = 0,
      int total_threads = -1);

  std::vector<std::vector<int>> GetMaxSimplices();
  friend void AddCofaces(const std::vector<std::vector<int>>& g, int depth,
                         int max_depth, std::vector<int> cur_node,
                         std::vector<int> neighbors, SimplicialComplex* simpl);

  /// binary[i][j] = 1 -> j-th simplex contains i element
  void BuildFromDowkerComplex(std::vector<std::vector<int>> binary);
};
