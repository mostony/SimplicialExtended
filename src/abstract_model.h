#pragma once

#include "Hasse.h"
#include <vector>
#include <string>
#include <functional>

class AbstractModel {
 public:
  void AddFunction(std::string name,
                   std::function<double(std::vector<int>)> func);
  void RemoveFunction(std::string name);
  std::vector<std::vector<double>> FeaturesMatrix(int rank);

  void ThresholdAbove(std::string name, double threshold);

  void ThresholdBelow(std::string name, double threshold);

  void UpdateWeight(std::vector<int> node, double new_weight);

  int HasseSize();

  std::vector<std::vector<int>> GetMaxFaces();
  std::vector<std::vector<int>> GetAll();
  std::vector<std::vector<int>> GetElementsWithRank(int rank);

  std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);
  int IncidenceDegree(std::vector<int> node, int k);

  std::vector<std::vector<int>> Adjacency(std::vector<int> node, int k);
  double Degree(std::vector<int> node, int k, bool weighted = false);
  std::vector<double> DegreeAll(int p, int k, bool weighted = false);

  int BettiNumber(int k);

  MyMatrixInt BoundaryMatrix(int k, int p);

  double Closeness(std::vector<int> node, int max_rank, bool weighted = false);
  double Betweenness(std::vector<int> node, int max_rank,
                     bool weighted = false);

  std::vector<std::pair<std::vector<int>, double>> ClosenessAll(
      int p, int max_rank, bool weighted = false);
  std::vector<std::pair<std::vector<int>, double>> BetweennessAll(
      int p, int max_rank, bool weighted = false);

  int Dimension();

  std::vector<std::pair<int, int>> FVector();

  int TotalCount();

  int EulerCharacteristic();

  void Clear();

  MyMatrixDiag Weights(int rank);

  MyMatrixDouble LaplacianMatrix(int k, int p, int q, bool weighted = false);
  std::vector<double> EigenValues(int k, int p, int q, bool weighted, int cnt);
  std::vector<double> EigenValuesAll(int k, int p, int q, bool weighted);

 protected:
  Hasse hasse_;
};
