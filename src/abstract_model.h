#pragma once

#include "Hasse.h"
#include <vector>

class AbstractModel {
public:
  int HasseSize();

  std::vector<std::vector<int>> GetMaxFaces();

  std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

  std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

  int BettiNumber(int k);

  double Closeness(std::vector<int> node, int max_rank);
  double Betweenness(std::vector<int> node, int max_rank);

protected:
  Hasse hasse_;
};
