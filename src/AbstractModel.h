#pragma once

#include <functional>
#include <string>
#include <vector>

#include "Hasse.h"

class AbstractModel {
   public:
    void AddFunction(std::string name, std::function<double(std::vector<int>)> func);
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

    std::vector<std::vector<int>> Adjacency(std::vector<int> node, int k);

    std::vector<std::vector<double>> AdjacencyMatrix(int k, int p, int q, bool weighted = false);
    double Degree(std::vector<int> node, int k, bool weighted = false);
    std::vector<double> DegreeAll(int p, int k, bool weighted = false);

    int BettiNumber(int k);

    double CommonNeighbors(std::vector<int> node1, std::vector<int> node2, int k);

    MyMatrixInt BoundaryMatrix(int k, int p);

    double Closeness(std::vector<int> node, int max_rank, bool weighted = false);
    double Betweenness(std::vector<int> node, int max_rank, bool weighted = false);

    std::vector<std::pair<std::vector<int>, double>> ClosenessAll(int p, int max_rank,
                                                                  bool weighted = false);
    std::vector<std::pair<std::vector<int>, double>> BetweennessAll(int p, int max_rank,
                                                                    bool weighted = false);

    std::vector<std::pair<std::vector<int>, double>> EigenCentrality(int p, int max_rank,
                                                                     bool weighted = false);
    std::vector<std::pair<std::vector<int>, double>> SubgraphCentrality(int p, int max_rank,
                                                                        bool weighted = false);
    int Dimension();

    std::vector<std::pair<int, int>> FVector();

    int TotalCount();

    int EulerCharacteristic();

    void Clear();

    MyMatrixDiag Weights(int rank);

    MyMatrixDouble LaplacianMatrix(int k, int p, int q, bool weighted = false,
                                   bool normalize = false);
    std::pair<std::vector<double>, std::vector<std::vector<double>>> EigenValues(
        int k, int p, int q, bool weighted, bool normalize, int cnt, const std::string& which);
    std::pair<std::vector<double>, std::vector<std::vector<double>>> EigenValuesAll(int k, int p,
                                                                                    int q,
                                                                                    bool weighted,
                                                                                    bool normalize);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> HodgeDecomposition(
        int k, int p, int q, const std::vector<double>& vec);

    /// on_column = true -- every column is simplex
    /// on_column = false -- every row is simplex
    virtual void BuildFromBinary(std::vector<std::vector<int>> binary, bool on_column) = 0;

   protected:
    Hasse hasse_;
};
