#include "abstract_model.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h"
#include "Hasse.h"

int AbstractModel::HasseSize() {
  return hasse_.Size();
}

std::vector<std::vector<int>> AbstractModel::GetMaxFaces() {
  return hasse_.GetMaxFaces();
}

std::vector<std::vector<int>> AbstractModel::GetAll() {
  return hasse_.GetAllElements();
}

std::vector<std::vector<int>> AbstractModel::Incidence(std::vector<int> node,
                                                       int k) {
  return hasse_.Incidence(node, k);
}

int AbstractModel::IncidenceDegree(std::vector<int> node, int k) {
  return hasse_.IncidenceDegree(node, k);
}

std::vector<std::vector<int>> AbstractModel::Adjacency(std::vector<int> node,
                                                       int k) {
  return hasse_.Adjacency(node, k);
}

int AbstractModel::Degree(std::vector<int> node, int k, bool weighted) {
  return hasse_.Degree(node, k, weighted);
}

int AbstractModel::BettiNumber(int k) {
  return hasse_.BettiNumber(k);
}

std::vector<std::vector<int>> AbstractModel::BoundaryMatrix(int k, int p) {
  auto result_mat = hasse_.BoundaryMatrix(k, p);
  size_t rows = result_mat.rows();
  size_t cols = result_mat.cols();
  std::vector<std::vector<int>> result(rows, std::vector<int>(cols));
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      result[i][j] = result_mat(i, j);
    }
  }
  return result;
}

double AbstractModel::Closeness(std::vector<int> node, int max_rank,
                                bool weighted) {
  return hasse_.Closeness(node, max_rank, weighted);
}

double AbstractModel::Betweenness(std::vector<int> node, int max_rank,
                                  bool weighted) {
  return hasse_.Betweenness(node, max_rank, weighted);
}

std::vector<std::pair<std::vector<int>, double>> AbstractModel::ClosenessAll(
    int p, int max_rank, bool weighted) {
  return hasse_.ClosenessAll(p, max_rank, weighted);
}

std::vector<std::pair<std::vector<int>, double>> AbstractModel::BetweennessAll(
    int p, int max_rank, bool weighted) {
  return hasse_.BetweennessAll(p, max_rank, weighted);
}

int AbstractModel::Dimension() {
  return hasse_.Dimension();
}

std::vector<std::pair<int, int>> AbstractModel::FVector() {
  return hasse_.FVector();
}

int AbstractModel::TotalCount() {
  return hasse_.TotalCount();
}

int AbstractModel::EulerCharacteristic() {
  return hasse_.EulerCharacteristic();
}

void AbstractModel::Clear() {
  hasse_ = Hasse();
}

std::vector<std::vector<double>> ConvertToVector(const MyMatrixDouble& mat) {
  std::vector<std::vector<double>> ret(mat.rows(),
                                       std::vector<double>(mat.cols()));
  for (size_t i = 0; i < mat.rows(); ++i) {
    for (size_t j = 0; j < mat.cols(); ++j) {
      ret[i][j] = mat(i, j);
    }
  }
  return ret;
}

std::vector<std::vector<double>> AbstractModel::Weights(int rank) {
  return ConvertToVector(hasse_.WeightedMatrix(rank));
}

std::vector<std::vector<double>> AbstractModel::LaplacianMatrix(int k, int p,
                                                                int q,
                                                                bool weighted) {
  auto L = hasse_.LaplacianMatrix(k, p, q, weighted);
  return ConvertToVector(L);
}

std::vector<double> AbstractModel::EigenValues(int k, int p, int q,
                                               bool weighted) {
  auto L = hasse_.LaplacianMatrix(k, p, q, weighted);
  Eigen::SelfAdjointEigenSolver<MyMatrixDouble> eg;
  auto values = eg.compute(L).eigenvalues();
  std::vector<double> result;
  for (auto x : values) {
    result.push_back(x);
  }
  return result;
}

void AbstractModel::AddFunction(std::string name,
                                std::function<double(std::vector<int>)> func) {
  hasse_.AddFunction(name, func);
}

void AbstractModel::RemoveFunction(std::string name) {
  return hasse_.RemoveFunction(name);
}

std::vector<std::vector<double>> AbstractModel::FeaturesMatrix(int rank) {
  return hasse_.FeaturesMatrix(rank);
}

void AbstractModel::ThresholdAbove(std::string name, double threshold) {
  hasse_.ThresholdAbove(name, threshold);
}

void AbstractModel::ThresholdBelow(std::string name, double threshold) {
  hasse_.ThresholdBelow(name, threshold);
}

void AbstractModel::UpdateWeight(std::vector<int> node, double new_weight) {
  hasse_.UpdateWeight(node, new_weight);
}
