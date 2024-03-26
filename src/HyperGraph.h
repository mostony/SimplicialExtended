#pragma once

#include <map>

#include "Hasse.h"

class HyperGraph {
public:
  void AddEdge(const std::vector<int> &edge);

  void RemoveEdge(const std::vector<int> &edge);

  std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

  std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

  int BettiNumber(int k);

  double Closeness(std::vector<int> node, int max_rank);
  double Betweenness(std::vector<int> node, int max_rank);

  std::vector<std::vector<int>> GetEdges();

private:
  Hasse hasse_;
};