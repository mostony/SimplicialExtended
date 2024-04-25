#pragma once

#include "Node.h"
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include <Eigen/Dense>

// long double throws bad::alloc on laplacian weighted for some reason
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MyMatrixDouble;

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MyMatrixInt;

class Hasse {
 public:
  /* ------------------- custom function ------------------ */
  void AddFunction(std::string name,
                   std::function<double(std::vector<int>)> func);

  void RemoveFunction(std::string name);

  std::vector<std::vector<double>> FeaturesMatrix(int rank);

  /* ----------------------- weights ---------------------- */
  void ThresholdAbove(std::string name, double threshold);

  void ThresholdBelow(std::string name, double threshold);

  void UpdateWeight(std::vector<int> node, double new_weight);
  // std::vector<std::vector<double>> GetWeights(int rank);

  /* ---------------------- structure --------------------- */
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

  /// return number of nodes
  int Size();

  /// Get max rank of structure
  int Dimension();

  /// merge two non intersecting Hasse diagrams
  friend void Merge(Hasse& current, Hasse& other);

  /* ----------------- Incidence & Degree ----------------- */

  /// only 0 and 1
  std::vector<std::vector<int>> IncidenceMatrix(int p, int k);

  // min weight
  std::vector<std::vector<double>> DegreeMatrix(int p, int k,
                                                bool weighted = false);

  std::vector<std::vector<int>> Incidence(const std::vector<int>& node, int k);
  int IncidenceDegree(const std::vector<int>& node, int k);

  std::vector<std::vector<int>> Adjacency(const std::vector<int>& node, int k);
  double Degree(const std::vector<int>& node, int k, bool weighted = false);

  int BettiNumber(int k);

  /* ---------------------- boundary ---------------------- */

  /// C_k -> C_p
  MyMatrixInt BoundaryMatrix(int k, int p);

  /// p < k < q
  MyMatrixDouble LaplacianMatrix(int k, int p, int q, bool weighted);

  /* --------------------- centrality --------------------- */
  double Closeness(std::vector<int> node, int max_rank, bool weighted = false);
  double Betweenness(std::vector<int> node, int max_rank,
                     bool weighted = false);

  /// get centrality stat for all nodes with rank = p
  /// through max_rank.
  std::vector<std::pair<std::vector<int>, double>> ClosenessAll(
      int p, int max_rank, bool weighted = false);
  std::vector<std::pair<std::vector<int>, double>> BetweennessAll(
      int p, int max_rank, bool weighted = false);

  std::vector<std::pair<int, int>> FVector();

  int TotalCount();

  int EulerCharacteristic();

  /// TODO: move this method out of class
  /// Check if a \in b
  static bool In(const std::vector<int>& a, const std::vector<int>& b);

  /// TODO: move this method out of class
  /// TODO: change int to enum
  /// O(n)
  static int CalculateSign(const std::vector<int>& subset,
                           const std::vector<int>& set);

  Eigen::DiagonalMatrix<double, Eigen::Dynamic> WeightedMatrix(int rank);

 private:
  void RecursiveRemoveNode(const std::vector<int>& node);

  void Threshold(std::function<bool(std::vector<int>)> is_good);

  void ResetCache();

  /// calculate SNF of mat module 2 and return
  /// number of non zero elements
  /// TODO: maybe some bitset optimizations?
  int SNF(std::vector<std::vector<int>> mat);

  std::vector<Node*> GetNodesWithFixedRank(int rank);

  int GetPositionInFixedRank(std::vector<int> node);

  std::map<std::vector<int>, std::unique_ptr<Node>> mapping_;

  std::map<int, std::set<Node*, decltype([](Node* a, Node* b) {
                           for (size_t i = 0;
                                i < std::min(a->data.size(), b->data.size());
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
  std::map<std::pair<int, int>, std::vector<std::vector<double>>> cache_degree_;
  std::map<std::string, std::function<double(std::vector<int>)>>
      custom_function_;
};
