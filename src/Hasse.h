#pragma once

#include "Node.h"
#include <map>
#include <memory>
#include <set>
#include <vector>

class Hasse {
public:
  void AddArc(const std::vector<int> &from, const std::vector<int> &to);

  /// remove node and all upper nodes, reachable from node
  void RemoveNode(const std::vector<int> &remove_node);

  void RemoveArc(const std::vector<int> &from, const std::vector<int> &to);

  /// add node and all lower nodes (subsets of node)
  void RecursiveAddNode(const std::vector<int> &add_node);

  /// use mapping to get node. If node not in mapping
  /// then create node. Return ptr to Node
  Node *GetNode(const std::vector<int> &node);

  std::vector<std::vector<int>> GetMaxFaces();

  std::vector<std::vector<int>> GetAllElements();

  /// return number of nodes
  int Size();

  /// merge two non intersecting Hasse diagrams
  friend void Merge(Hasse &current, Hasse &other);

  std::vector<std::vector<int>> IncidenceMatrix(int p, int k);
  std::vector<std::vector<int>> DegreeMatrix(int p, int k);

  /// calculate SNF of mat module 2 and return
  /// number of non zero elements
  /// TODO: maybe some bitset optimizations?
  int SNF(std::vector<std::vector<int>> mat);

  int BettiNumber(int k);

  std::vector<std::vector<int>> Incidence(const std::vector<int> &node, int k);
  int IncidenceDegree(const std::vector<int> &node, int k);

  /// TODO: maybe remove self node from degree???
  std::vector<std::vector<int>> Adjacency(const std::vector<int> &node, int k);
  int Degree(const std::vector<int> &node, int k);

  double Closeness(std::vector<int> node, int max_rank);
  double Betweenness(std::vector<int> node, int max_rank);

  /// Check if a \in b
  static bool In(const std::vector<int> &a, const std::vector<int> &b);

private:
  void RecursiveRemoveNode(const std::vector<int> &node);

  std::map<std::vector<int>, std::unique_ptr<Node>> mapping_;

  std::map<int, std::set<Node *>> nodes_with_fixed_rank_;

  /// TODO: fix when removing node
  /// in this case can't use vector or need to rebuild
  /// 1) add flag to cache
  /// 2) change nodes_with_fixed_rank_ (when no edges to it)
  std::map<std::pair<int, int>, std::vector<std::vector<int>>> cache_incidence_;
  std::map<std::pair<int, int>, std::vector<std::vector<int>>> cache_degree_;
};
