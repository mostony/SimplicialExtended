#include "CombinatorialComplex.h"

int CombinatorialComplex::HasseSize() { return hasse_.Size(); }

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

std::vector<std::vector<int>>
CombinatorialComplex::Incidence(std::vector<int> node, int k) {
  return hasse_.Incidence(node, k);
}

std::vector<std::vector<int>>
CombinatorialComplex::Degree(std::vector<int> node, int k) {
  return hasse_.Degree(node, k);
}

int CombinatorialComplex::BettiNumber(int k) { return hasse_.BettiNumber(k); }

double CombinatorialComplex::Closeness(std::vector<int> node, int max_rank) {
  return hasse_.Closeness(node, max_rank);
}

double CombinatorialComplex::Betweenness(std::vector<int> node, int max_rank) {
  return hasse_.Betweenness(node, max_rank);
}

std::vector<std::vector<int>> CombinatorialComplex::GetSubsets() {
  return hasse_.GetMaxFaces();
}
