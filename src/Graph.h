#pragma once

#include "Hasse.h"

class Graph {
public:
  void AddEdge(int v, int u);

  void RemoveEdge(int v, int u);

  std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

  std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

  int BettiNumber(int k);

  double Closeness(std::vector<int> node, int max_rank);

  double Betweenness(std::vector<int> node, int max_rank);

private:
  Hasse hasse_;
};