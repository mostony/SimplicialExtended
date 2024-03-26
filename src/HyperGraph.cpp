#include "HyperGraph.h"

void HyperGraph::AddEdge(const std::vector<int> &edge) {
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.AddArc(lower_node, upper_node);
  }
}

void HyperGraph::RemoveEdge(const std::vector<int> &edge) {
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.RemoveArc(lower_node, upper_node);
  }
}

std::vector<std::vector<int>> HyperGraph::Incidence(std::vector<int> node,
                                                    int k) {
  return hasse_.Incidence(node, k);
}

std::vector<std::vector<int>> HyperGraph::Degree(std::vector<int> node, int k) {
  return hasse_.Degree(node, k);
}

int HyperGraph::BettiNumber(int k) { return hasse_.BettiNumber(k); }

double HyperGraph::Closeness(std::vector<int> node, int max_rank) {
  return hasse_.Closeness(node, max_rank);
}

double HyperGraph::Betweenness(std::vector<int> node, int max_rank) {
  return hasse_.Closeness(node, max_rank);
}

std::vector<std::vector<int>> HyperGraph::GetEdges() {
  return hasse_.GetMaxFaces();
}
