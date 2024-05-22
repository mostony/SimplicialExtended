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

void HyperGraph::BuildFromBinary(std::vector<std::vector<int>> binary,
                                 bool on_column) {
  Clear();
  if (binary.empty()) {
    return;
  }

  if (on_column == false) {
    for (size_t i = 0; i < binary.size(); i++) {
      std::vector<int> edge;
      for (size_t j = 0; j < binary[i].size(); j++) {
        if (binary[i][j]) {
          edge.push_back(j);
        }
      }
      AddEdge(edge);
    }
  } else {
    for (size_t j = 0; j < binary[0].size(); j++) {
      std::vector<int> edge;
      for (size_t i = 0; i < binary.size(); i++) {
        if (binary[i][j]) {
          edge.push_back(i);
        }
      }
      AddEdge(edge);
    }
  }
}
