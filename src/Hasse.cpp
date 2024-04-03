#include "Hasse.h"

#include <cassert>
#include <queue>
#include <set>
#include <stdexcept>

Node *Hasse::GetNode(const std::vector<int> &node) {
  if (!mapping_.count(node)) {
    mapping_[node] = std::unique_ptr<Node>(new Node(node));
    nodes_with_fixed_rank_[mapping_[node]->rank].insert(mapping_[node].get());
  }
  return mapping_[node].get();
}

void Hasse::AddArc(const std::vector<int> &from, const std::vector<int> &to) {
  auto low = GetNode(from);
  auto up = GetNode(to);
  low->upper.push_back(to);
  up->lower.push_back(from);
}

/// TODO: add remove from mapping
void Hasse::RemoveNode(const std::vector<int> &node) {
  RecursiveRemoveNode(node);
}

void Hasse::RemoveArc(const std::vector<int> &from,
                      const std::vector<int> &to) {
  auto low = GetNode(from);
  auto up = GetNode(to);
  for (auto it = low->upper.begin(); it != low->upper.end(); it++) {
    if (*it == to) {
      low->upper.erase(it);
      break;
    }
  }

  for (auto it = up->lower.begin(); it != up->lower.end(); it++) {
    if (*it == from) {
      up->lower.erase(it);
      break;
    }
  }
}

void Hasse::RecursiveRemoveNode(const std::vector<int> &remove_node) {
  std::queue<std::vector<int>> q;
  std::set<std::vector<int>> used;
  q.push(remove_node);
  used.insert(remove_node);
  while (!q.empty()) {
    auto top = q.front();
    q.pop();
    auto node = GetNode(top);
    nodes_with_fixed_rank_[node->rank].erase(node);

    for (const auto &nxt : node->upper) {
      if (!used.count(nxt)) {
        used.insert(nxt);
        q.push(nxt);
      }
    }

    while (!node->lower.empty()) {
      const auto &prev = node->lower.back();
      RemoveArc(prev, top);
    }

    node->upper.clear();
    node->lower.clear();
  }

  for (const auto &node : used) {
    mapping_.erase(node);
  }
}

bool Hasse::In(const std::vector<int> &a, const std::vector<int> &b) {
  for (size_t i = 0; i < a.size(); i++) {
    bool found = false;
    for (size_t j = 0; j < b.size(); j++) {
      if (a[i] == b[j]) {
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

void Hasse::RecursiveAddNode(const std::vector<int> &add_node) {
  if (mapping_.count(add_node)) {
    return;
  }

  std::queue<std::vector<int>> q;
  q.push(add_node);
  while (!q.empty()) {
    auto top = q.front();
    q.pop();

    auto node = GetNode(top);
    if (node->size == 0) {
      continue;
    }

    for (size_t i = 0; i < node->size; i++) {
      auto next = top;
      next.erase(next.begin() + i);

      bool was_before = mapping_.count(next);
      auto low = GetNode(next);

      AddArc(low->data, top);

      if (!was_before) {
        q.push(next);
      }
    }
  }
}

std::vector<std::vector<int>> Hasse::GetMaxFaces() {
  std::vector<std::vector<int>> result;
  for (auto &[id, node] : mapping_) {
    if (node->upper.empty()) {
      result.push_back(node->data);
    }
  }
  return result;
}

std::vector<std::vector<int>> Hasse::GetAllElements() {
  std::vector<std::vector<int>> result;
  for (auto &[id, node] : mapping_) {
    result.push_back(node->data);
  }
  return result;
}

int Hasse::Size() { return mapping_.size(); }

std::vector<std::vector<int>> Hasse::IncidenceMatrix(int p, int k) {
  if (k < p) {
    throw std::runtime_error("k < size of node");
  }
  // p->k
  if (cache_incidence_.count(std::make_pair(k, p))) {
    return cache_incidence_[std::make_pair(k, p)];
  }
  std::vector<Node *> lower(nodes_with_fixed_rank_[p].begin(),
                            nodes_with_fixed_rank_[p].end());
  std::vector<Node *> upper(nodes_with_fixed_rank_[k].begin(),
                            nodes_with_fixed_rank_[k].end());
  std::vector<std::vector<int>> result(lower.size(),
                                       std::vector<int>(upper.size()));
  for (size_t i = 0; i < lower.size(); i++) {
    for (size_t j = 0; j < upper.size(); j++) {
      result[i][j] = (Hasse::In(lower[i]->data, upper[j]->data));
    }
  }
  cache_incidence_[std::make_pair(k, p)] = std::move(result);
  return cache_incidence_[std::make_pair(k, p)];
}

std::vector<std::vector<int>> Hasse::DegreeMatrix(int p, int k) {
  if (k < p) {
    throw std::runtime_error("k < size of node");
  }
  if (cache_degree_.count(std::make_pair(k, p))) {
    return cache_degree_[std::make_pair(k, p)];
  }

  auto incidence = IncidenceMatrix(p, k);
  if (incidence.empty()) {
    return {};
  }

  int size = incidence.size();
  int inter_size = incidence[0].size();
  std::vector<std::vector<int>> result(size, std::vector<int>(size));
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      for (size_t k = 0; k < inter_size; k++) {
        if (incidence[i][k] && incidence[j][k]) {
          result[i][j] = true;
          break;
        }
      }
    }
  }
  cache_degree_[std::make_pair(k, p)] = std::move(result);
  return cache_degree_[std::make_pair(k, p)];
}

int Hasse::SNF(std::vector<std::vector<int>> mat) {
  if (mat.empty()) {
    return 0;
  }
  for (auto row : mat) {
    for (auto elem : row) {
      // modulo 2
      assert(0 <= elem && elem <= 1);
    }
  }
  int rows = mat.size();
  int cols = mat[0].size();
  for (size_t k = 0; k < std::min(rows, cols); k++) {
    bool found = false;
    for (size_t i = k; i < rows; i++) {
      if (found) {
        break;
      }
      for (size_t j = k; j < cols; j++) {
        if (!mat[i][j]) {
          continue;
        }

        // swap i and k rows
        if (i != k) {
          for (size_t y = k; y < cols; y++) {
            std::swap(mat[i][y], mat[k][y]);
          }
        }

        // swap j and k columns
        if (j != k) {
          for (size_t x = k; x < rows; x++) {
            std::swap(mat[x][j], mat[x][k]);
          }
        }

        // add row
        for (size_t x = k + 1; x < rows; x++) {
          if (mat[x][k]) {
            for (size_t y = k; y < cols; y++) {
              mat[x][y] ^= mat[k][y];
            }
          }
        }

        // add column
        for (size_t y = k + 1; y < cols; y++) {
          if (mat[k][y]) {
            for (size_t x = k; x < rows; x++) {
              mat[x][y] ^= mat[x][k];
            }
          }
        }
        found = true;
        break;
      }
    }
    if (!found) {
      return k;
    }
  }
  return std::min(rows, cols);
}

int Hasse::BettiNumber(int k) {
  /// TODO: add runtime checks
  /// TODO: check code
  int kernel = 0;
  if (k == 0) {
    kernel = this->nodes_with_fixed_rank_[0].size();
  } else {
    auto inc = IncidenceMatrix(k - 1, k);
    if (!inc.empty()) {
      kernel = (int)(inc[0].size()) - SNF(inc);
    }
  }
  int image = 0;
  if (nodes_with_fixed_rank_[k + 1].size() != 0) {
    auto inc = IncidenceMatrix(k, k + 1);
    image = SNF(inc);
  }
  return kernel - image;
}

void Merge(Hasse &current, Hasse &other) {
  for (auto &[key, value] : other.mapping_) {
    current.mapping_[key] = std::move(value);
  }
}

std::vector<std::vector<int>> Hasse::Incidence(const std::vector<int> &node,
                                               int k) {
  int p = GetNode(node)->rank;
  if (k < p) {
    throw std::runtime_error("k < size of node");
  }
  auto mat = IncidenceMatrix(p, k);
  std::vector<std::vector<int>> result;
  std::vector<Node *> upper(nodes_with_fixed_rank_[k].begin(),
                            nodes_with_fixed_rank_[k].end());
  std::vector<Node *> lower(nodes_with_fixed_rank_[p].begin(),
                            nodes_with_fixed_rank_[p].end());
  size_t row = 0;
  while (lower[row]->data != node) {
    row += 1;
  }
  // TODO: make bitset or just links
  for (size_t col = 0; col < mat[row].size(); col++) {
    if (mat[row][col]) {
      result.push_back(upper[col]->data);
    }
  }
  return result;
}

std::vector<std::vector<int>> Hasse::Degree(const std::vector<int> &node,
                                            int k) {
  int p = GetNode(node)->rank;
  std::vector<std::vector<int>> result;
  auto mat = DegreeMatrix(p, k);

  std::vector<Node *> simplices(nodes_with_fixed_rank_[p].begin(),
                                nodes_with_fixed_rank_[p].end());
  size_t row = 0;
  while (simplices[row]->data != node) {
    row += 1;
  }

  for (size_t col = 0; col < mat[row].size(); col++) {
    if (mat[row][col]) {
      result.push_back(simplices[col]->data);
    }
  }
  return result;
}

double Hasse::Closeness(std::vector<int> node, int max_rank) {
  if (!mapping_.count(node)) {
    throw std::runtime_error("No such node");
  }
  auto f = GetNode(node);
  int k = f->rank;
  if (max_rank <= k) {
    throw std::runtime_error("k >= max_rank");
  }
  if (nodes_with_fixed_rank_[k].size() == 1) {
    throw std::runtime_error("no other nodes");
  }

  auto g = DegreeMatrix(k, max_rank);
  int n = g.size();

  std::vector<int> dist(n, n);
  std::queue<int> q;

  // TODO: very bad code: add function
  std::vector<Node *> tmp(nodes_with_fixed_rank_[k].begin(),
                          nodes_with_fixed_rank_[k].end());
  size_t start = 0;
  while (tmp[start]->data != node) {
    start += 1;
  }
  q.push(start);
  dist[start] = 0;
  while (!q.empty()) {
    int v = q.front();
    q.pop();
    for (int u = 0; u < n; u++) {
      if (g[v][u] && dist[u] > dist[v] + 1) {
        dist[u] = dist[v] + 1;
        q.push(u);
      }
    }
  }
  int sum_distances = 0;
  bool is_connected = true;
  for (int v = 0; v < n; v++) {
    sum_distances += dist[v];
    if (v != start && dist[v] == n) {
      is_connected = false;
    }
  }
  // TODO: if not connected equal zero ???
  if (!is_connected) {
    return 0.0;
  }
  double avg_sum_distance = 1.0 * sum_distances / (n - 1);
  return 1 / avg_sum_distance;
}

double Hasse::Betweenness(std::vector<int> node, int max_rank) {
  if (!mapping_.count(node)) {
    throw std::runtime_error("No such node");
  }
  auto f = GetNode(node);
  int k = f->rank;
  if (max_rank <= k) {
    throw std::runtime_error("k >= max_rank");
  }
  if (nodes_with_fixed_rank_[k].size() <= 2) {
    throw std::runtime_error("too low other nodes");
  }

  auto g = DegreeMatrix(k, max_rank);
  int n = g.size();

  std::vector<Node *> nodes(nodes_with_fixed_rank_[k].begin(),
                            nodes_with_fixed_rank_[k].end());
  size_t ind_node = 0;
  while (nodes[ind_node]->data != node) {
    ind_node += 1;
  }
  double sum_distances = 0;

  for (size_t i = 0; i < nodes.size(); i++) {
    auto s = nodes[i];
    if (i == ind_node) {
      continue;
    }
    std::vector<int> dist(n, n);
    std::vector<int> shortest_path_cnt(n, 0);
    std::vector<std::pair<int, int>> shortest_path_through_node_cnt(
        n, {0, 0}); // (cnt, visited)
    std::queue<int> q;
    dist[i] = 0;
    q.push(i);
    shortest_path_cnt[i] = 1;
    shortest_path_through_node_cnt[i] = {1, 0};

    while (!q.empty()) {
      int v = q.front();
      q.pop();
      if (v == ind_node) {
        shortest_path_through_node_cnt[v].second = true;
      }
      for (int u = 0; u < n; u++) {
        if (!g[v][u]) {
          continue;
        }
        if (dist[u] > dist[v] + 1) {
          dist[u] = dist[v] + 1;
          shortest_path_cnt[u] = shortest_path_cnt[v];
          shortest_path_through_node_cnt[u] = shortest_path_through_node_cnt[v];
          q.push(u);
        } else if (dist[u] == dist[v] + 1) {
          shortest_path_cnt[u] += shortest_path_cnt[v];

          if (shortest_path_through_node_cnt[u].second == 0 &&
              shortest_path_through_node_cnt[v].second == 1) {
            shortest_path_through_node_cnt[u] =
                shortest_path_through_node_cnt[v];
          } else if (shortest_path_through_node_cnt[u].second == 1 &&
                     shortest_path_through_node_cnt[v].second == 0) {
            // skip
          } else {
            shortest_path_through_node_cnt[u].first +=
                shortest_path_through_node_cnt[v].first;
          }
        }
      }
    }
    for (int u = i + 1; u < n; u++) {
      if (u == ind_node) {
        continue;
      }
      int num = shortest_path_through_node_cnt[u].first *
                shortest_path_through_node_cnt[u].second;
      int den = shortest_path_cnt[u];
      if (den == 0) {
        throw std::runtime_error("Unconnected graph");
      }
      sum_distances += 1.0 * num / den;
    }
  }
  int norm = (n - 1) * (n - 2) / 2;
  return sum_distances / norm;
}