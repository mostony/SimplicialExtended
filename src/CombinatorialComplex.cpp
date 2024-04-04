#include "CombinatorialComplex.h"

void CombinatorialComplex::Build(const std::vector<std::vector<int>> &data) {
  std::vector<Node *> nodes;
  for (size_t i = 0; i < data.size(); i++) {
    nodes.push_back(hasse_.GetNode(data[i]));
  }

  for (size_t i = 0; i < data.size(); i++) {
    for (size_t j = 0; j < data.size(); j++) {
      if (i == j) {
        continue;
      }
      if (Hasse::In(data[i], data[j])) {
        hasse_.AddArc(data[i], data[j]);
      }
    }
  }
}

std::vector<std::vector<int>> CombinatorialComplex::GetSubsets() {
  return GetMaxFaces();
}
