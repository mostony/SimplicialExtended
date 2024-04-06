#include "abstract_model.h"
#include "Hasse.h"

int AbstractModel::HasseSize() {
  return hasse_.Size();
}

std::vector<std::vector<int>> AbstractModel::GetMaxFaces() {
  return hasse_.GetMaxFaces();
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

int AbstractModel::Degree(std::vector<int> node, int k) {
  return hasse_.Degree(node, k);
}

int AbstractModel::BettiNumber(int k) {
  return hasse_.BettiNumber(k);
}

double AbstractModel::Closeness(std::vector<int> node, int max_rank) {
  return hasse_.Closeness(node, max_rank);
}

double AbstractModel::Betweenness(std::vector<int> node, int max_rank) {
  return hasse_.Betweenness(node, max_rank);
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
