#include "Graph.h"
#include "Hasse.h"

void Graph::AddEdge(int v, int u) {
    hasse_.RecursiveAddNode({v, u});
}

void Graph::RemoveEdge(int v, int u) {
    hasse_.RemoveNode({v, u});
}

std::vector<std::vector<int>> Graph::Incidence(std::vector<int> node, int k) {
    return hasse_.Incidence(node, k);
}

std::vector<std::vector<int>> Graph::Degree(std::vector<int> node, int k) {
    return hasse_.Degree(node, k);
}

int Graph::BettiNumber(int k) {
    return hasse_.BettiNumber(k);
}

double Graph::Closeness(std::vector<int> node, int max_rank) {
    return hasse_.Closeness(node, max_rank);
}
