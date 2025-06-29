#pragma once

#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include "Node.h"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> MyMatrixDouble;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> MyMatrixDiag;

typedef Eigen::SparseMatrix<int, Eigen::RowMajor> MyMatrixInt;

class Hasse {
   public:
    /* ------------------- custom function ------------------ */
    void AddFunction(std::string name, std::function<double(std::vector<int>)> func);

    void RemoveFunction(std::string name);

    std::vector<std::vector<double>> FeaturesMatrix(int rank);

    /* ----------------------- weights ---------------------- */
    void ThresholdAbove(std::string name, double threshold);

    void ThresholdBelow(std::string name, double threshold);

    void UpdateWeight(std::vector<int> node, double new_weight);

    /* ---------------------- structure --------------------- */

    void CreateNode(const std::vector<int>& data, int rank);

    void AddArc(const std::vector<int>& from, const std::vector<int>& to);

    /// remove node and all upper nodes, reachable from node
    void RemoveNode(const std::vector<int>& remove_node);

    void RemoveArc(const std::vector<int>& from, const std::vector<int>& to);

    /// add node and all subsets of node
    void RecursiveAddNode(const std::vector<int>& add_node);

    Node* GetNode(const std::vector<int>& node);

    /* ------------------- representation ------------------- */
    std::vector<std::vector<int>> GetMaxFaces();

    std::vector<std::vector<int>> GetAllElements();

    std::vector<std::vector<int>> GetElementsWithRank(int rank);

    /// return number of nodes
    int Size();

    /// Get max rank of structure
    int Dimension();

    /// merge two non intersecting Hasse diagrams
    friend void Merge(Hasse& current, Hasse& other);

    /* ----------------- Incidence & Degree ----------------- */

    /// only 0 and 1
    std::vector<std::vector<int>>& IncidenceMatrix(int p, int k);

    /// Min weight
    std::vector<std::vector<double>>& UpperAdjacencyMatrix(int p, int k, bool weighted = false);

    /// Min weight
    std::vector<std::vector<double>> LowerAdjacencyMatrix(int p, int k, bool weighted = false);

    /// Adjacency matrix equals difference between lower adjacency and upper adjacency matrixes.
    /// In other words: (i, j) == 1 iff i_th and j_th simplices of rank k are lower adjacent, but
    /// not upper adjacent. But k==0 is special case in which i_th and j_th simplices if they upper
    /// adjacent.
    std::vector<std::vector<double>> AdjacencyMatrix(int k, int p, int q, bool weighted = false);

    std::vector<std::vector<int>> Incidence(const std::vector<int>& node, int k);

    std::vector<std::vector<int>> Adjacency(const std::vector<int>& node, int k);
    double Degree(const std::vector<int>& node, int k, bool weighted = false);

    std::vector<double> DegreeAll(int p, int k, bool weighted = false);

    int BettiNumber(int k);

    double CommonNeighbors(std::vector<int> node1, std::vector<int> node2, int k);

    /* ---------------------- boundary ---------------------- */

    /// C_k -> C_p
    MyMatrixInt BoundaryMatrix(int k, int p);

    /// p < k < q
    MyMatrixDouble LaplacianMatrix(int k, int p, int q, bool weighted, bool normalize);

    std::pair<std::vector<double>, std::vector<std::vector<double>>> EigenValues(
        int k, int p, int q, bool weighted, bool normalize, int cnt, const std::string& which);
    std::pair<std::vector<double>, std::vector<std::vector<double>>> EigenValuesAll(int k, int p,
                                                                                    int q,
                                                                                    bool weighted,
                                                                                    bool normalize);

    /// Returns gradient, harmonic and solenoidal components
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> HodgeDecomposition(
        int k, int p, int q, const std::vector<double>& vec);

    /* --------------------- centrality --------------------- */
    double Closeness(std::vector<int> node, int max_rank, bool weighted = false);
    double Betweenness(std::vector<int> node, int max_rank, bool weighted = false);

    /// get centrality stat for all nodes with rank = p
    /// through max_rank.
    std::vector<std::pair<std::vector<int>, double>> ClosenessAll(int p, int max_rank,
                                                                  bool weighted = false);
    std::vector<std::pair<std::vector<int>, double>> BetweennessAll(int p, int max_rank,
                                                                    bool weighted = false);

    std::vector<std::pair<std::vector<int>, double>> EigenCentrality(int p, int max_rank,
                                                                     bool weighted = false);
    std::vector<std::pair<std::vector<int>, double>> SubgraphCentrality(int p, int max_rank,
                                                                        bool weighted = false);

    std::vector<std::pair<int, int>> FVector();

    int TotalCount();

    int EulerCharacteristic();

    static int CalculateSign(const std::vector<int>& subset, const std::vector<int>& set);

    MyMatrixDiag WeightedMatrix(int rank);

    MyMatrixDiag DegreeMatrix(int p, int k, bool weighted = false);

   private:
    void RecursiveRemoveNode(const std::vector<int>& node);

    void Threshold(std::function<bool(std::vector<int>)> is_good);

    void ResetCache();

    /// calculate SNF of mat module 2 and return
    /// number of non zero elements
    int SNF(std::vector<std::vector<int>> mat);

    std::vector<Node*> GetNodesWithFixedRank(int rank);

    int GetPositionInFixedRank(std::vector<int> node);

    // Projection on columns
    static Eigen::VectorXd OrthProject(const std::vector<double>& vec, const MyMatrixDouble& A);

    std::map<std::vector<int>, std::unique_ptr<Node>> mapping_;

    std::map<int, std::set<Node*, decltype([](Node* a, Node* b) {
                               for (size_t i = 0; i < std::min(a->data.size(), b->data.size());
                                    i++) {
                                   if (a->data[i] < b->data[i]) {
                                       return true;
                                   }
                                   if (a->data[i] > b->data[i]) {
                                       return false;
                                   }
                               }
                               return a->data.size() < b->data.size();
                           })>>
        nodes_with_fixed_rank_;

    std::map<std::pair<int, int>, std::vector<std::vector<int>>> cache_incidence_;
    std::map<std::array<int, 3>, std::vector<std::vector<double>>> cache_degree_;
    std::map<std::pair<int, int>, MyMatrixInt> cache_boundary_;
    std::map<std::string, std::function<double(std::vector<int>)>> custom_function_;

    static constexpr double EPS = 1e-9;
};

bool In(const std::vector<int>& a, const std::vector<int>& b);
