#pragma once

#include "Hasse.h"
#include <vector>

class AbstractModel {
 public:
  void AddFunction(std::string name,
                   std::function<double(std::vector<int>)> func);
  void ThresholdAbove(std::string name, double threshold);

  int HasseSize();

  std::vector<std::vector<int>> GetMaxFaces();

  std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);
  int IncidenceDegree(std::vector<int> node, int k);

  std::vector<std::vector<int>> Adjacency(std::vector<int> node, int k);
  int Degree(std::vector<int> node, int k);

  int BettiNumber(int k);

  double Closeness(std::vector<int> node, int max_rank);
  double Betweenness(std::vector<int> node, int max_rank);

  std::vector<std::pair<std::vector<int>, double>> ClosenessAll(int p,
                                                                int max_rank);
  std::vector<std::pair<std::vector<int>, double>> BetweennessAll(int p,
                                                                  int max_rank);

  int Dimension();

  std::vector<std::pair<int, int>> FVector();

  int TotalCount();

  int EulerCharacteristic();

  void Clear();

 protected:
  Hasse hasse_;
};
