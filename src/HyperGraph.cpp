#include "HyperGraph.h"

void HyperGraph::AddEdge(const std::vector<int>& edge) {
    for (auto vertice : edge) {
        auto lower_node = std::vector<int>{vertice};
        auto upper_node = edge;
        hasse_.AddArc(lower_node, upper_node);
    }
}

void HyperGraph::RemoveEdge(const std::vector<int>& edge) {
    for (auto vertice : edge) {
        auto lower_node = std::vector<int>{vertice};
        auto upper_node = edge;
        hasse_.RemoveArc(lower_node, upper_node);
    }
}
