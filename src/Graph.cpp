#include "Graph.h"
#include "Hasse.h"

void Graph::AddEdge(int v, int u) {
  if (v > u) {
    std::swap(v, u);
  }
  hasse_.RecursiveAddNode({v, u});
}

void Graph::RemoveEdge(int v, int u) {
  if (v > u) {
    std::swap(v, u);
  }
  hasse_.RemoveNode({v, u});
}

std::vector<std::vector<int>> Graph::GetEdges() {
  return hasse_.GetMaxFaces();
}

void Graph::BuildFromBinary(std::vector<std::vector<int>> binary,
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
      for (auto v : edge) {
        for (auto u : edge) {
          if (v < u) {
            AddEdge(v, u);
          }
        }
      }
    }
  } else {
    for (size_t j = 0; j < binary[0].size(); j++) {
      std::vector<int> edge;
      for (size_t i = 0; i < binary.size(); i++) {
        if (binary[i][j]) {
          edge.push_back(i);
        }
      }
      for (auto v : edge) {
        for (auto u : edge) {
          if (v < u) {
            AddEdge(v, u);
          }
        }
      }
    }
  }
}
