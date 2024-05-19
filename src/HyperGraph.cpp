#include "HyperGraph.h"
#include <algorithm>

void HyperGraph::AddEdge(std::vector<int> edge) {
  std::sort(edge.begin(), edge.end());
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.CreateNode(upper_node, 1);
    hasse_.CreateNode(lower_node, 0);

    hasse_.AddArc(lower_node, upper_node);
  }
}

void HyperGraph::RemoveEdge(std::vector<int> edge) {
  std::sort(edge.begin(), edge.end());
  for (auto vertice : edge) {
    auto lower_node = std::vector<int>{vertice};
    auto upper_node = edge;
    hasse_.RemoveArc(lower_node, upper_node);
  }
}

std::vector<std::vector<int>> HyperGraph::GetEdges() {
  return this->GetMaxFaces();
}
