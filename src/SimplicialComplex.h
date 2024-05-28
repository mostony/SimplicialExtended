#pragma once

#include <vector>

#include "AbstractModel.h"

class SimplicialComplex : public AbstractModel {
 public:
  void AddSimplex(std::vector<int> simplex);

  void RemoveSimplex(std::vector<int> simplex);

  /// merge two non intersecting SimplicialComplex
  friend void Merge(SimplicialComplex* current, SimplicialComplex* other);

  /// TODO: change method to enum
  static std::unique_ptr<SimplicialComplex> CreateCliqueGraph(
      const std::vector<std::vector<int>>& g, int k, int method = 0,
      int total_threads = -1);

  std::vector<std::vector<int>> GetMaxSimplices();
  friend void AddCofaces(const std::vector<std::vector<int>>& g, int depth,
                         int max_depth, std::vector<int> cur_node,
                         std::vector<int> neighbors, SimplicialComplex* simpl);

  void BuildFromBinary(std::vector<std::vector<int>> binary,
                       bool on_column) override;
};
