#include "HyperGraph.h"

void HyperGraph::AddEdge(std::vector<int> edge) {
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.AddArc(lower_node, upper_node);
  }
}

void HyperGraph::RemoveEdge(std::vector<int> edge) {
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.RemoveArc(lower_node, upper_node);
  }
}

std::vector<std::vector<int>> HyperGraph::GetEdges() {
  return this->GetMaxFaces();
}
