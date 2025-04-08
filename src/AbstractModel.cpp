#include "AbstractModel.h"

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

std::vector<std::vector<int>> AbstractModel::Incidence(std::vector<int> node, int k) {
    return hasse_.Incidence(node, k);
}

int AbstractModel::IncidenceDegree(std::vector<int> node, int k) {
    return hasse_.IncidenceDegree(node, k);
}

std::vector<std::vector<int>> AbstractModel::Adjacency(std::vector<int> node, int k) {
    return hasse_.Adjacency(node, k);
}

double AbstractModel::Degree(std::vector<int> node, int k, bool weighted) {
    return hasse_.Degree(node, k, weighted);
}

std::vector<double> AbstractModel::DegreeAll(int p, int k, bool weighted) {
    return hasse_.DegreeAll(p, k, weighted);
}

int AbstractModel::BettiNumber(int k) {
    return hasse_.BettiNumber(k);
}

MyMatrixInt AbstractModel::BoundaryMatrix(int k, int p) {
    return hasse_.BoundaryMatrix(k, p);
}

double AbstractModel::Closeness(std::vector<int> node, int max_rank, bool weighted) {
    return hasse_.Closeness(node, max_rank, weighted);
}

double AbstractModel::Betweenness(std::vector<int> node, int max_rank, bool weighted) {
    return hasse_.Betweenness(node, max_rank, weighted);
}

std::vector<std::pair<std::vector<int>, double>> AbstractModel::ClosenessAll(int p, int max_rank,
                                                                             bool weighted) {
    return hasse_.ClosenessAll(p, max_rank, weighted);
}

std::vector<std::pair<std::vector<int>, double>> AbstractModel::BetweennessAll(int p, int max_rank,
                                                                               bool weighted) {
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
    std::vector<std::vector<double>> ret(mat.rows(), std::vector<double>(mat.cols()));
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            ret[i][j] = mat.coeff(i, j);
        }
    }
    return ret;
}

std::vector<std::vector<double>> ConvertToVectorDiag(const MyMatrixDiag& mat) {
    std::vector<std::vector<double>> ret(mat.rows(), std::vector<double>(mat.cols()));
    for (size_t i = 0; i < mat.rows(); ++i) {
        ret[i][i] = mat.diagonal()[i];
    }
    return ret;
}

MyMatrixDiag AbstractModel::Weights(int rank) {
    return hasse_.WeightedMatrix(rank);
}

MyMatrixDouble AbstractModel::LaplacianMatrix(int k, int p, int q, bool weighted) {
    return hasse_.LaplacianMatrix(k, p, q, weighted);
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> AbstractModel::EigenValues(
    int k, int p, int q, bool weighted, int cnt, const std::string& which) {
    return hasse_.EigenValues(k, p, q, weighted, cnt, which);
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> AbstractModel::EigenValuesAll(
    int k, int p, int q, bool weighted) {
    return hasse_.EigenValuesAll(k, p, q, weighted);
}

void AbstractModel::AddFunction(std::string name, std::function<double(std::vector<int>)> func) {
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

std::vector<std::vector<int>> AbstractModel::GetElementsWithRank(int rank) {
    return hasse_.GetElementsWithRank(rank);
}
