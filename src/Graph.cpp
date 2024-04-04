#include "Graph.h"
#include "Hasse.h"

void Graph::AddEdge(int v, int u) { hasse_.RecursiveAddNode({v, u}); }

void Graph::RemoveEdge(int v, int u) { hasse_.RemoveNode({v, u}); }

std::vector<std::vector<int>> Graph::GetEdges() { return hasse_.GetMaxFaces(); }
