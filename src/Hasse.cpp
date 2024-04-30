#include "Hasse.h"
#include "Node.h"

#include <cassert>
#include <functional>
#include <mutex>
#include <queue>
#include <set>
#include <stdexcept>

#include <algorithm>
#include <thread>
#include <unordered_map>
#include <vector>

Node* Hasse::GetNode(const std::vector<int>& node) {
  assert(is_sorted(node.begin(), node.end()));
  if (!mapping_.count(node)) {
    mapping_[node] = std::unique_ptr<Node>(new Node(node));
    nodes_with_fixed_rank_[mapping_[node]->rank].insert(mapping_[node].get());
  }
  return mapping_[node].get();
}

void Hasse::AddArc(const std::vector<int>& from, const std::vector<int>& to) {
  ResetCache();
  auto low = GetNode(from);
  auto up = GetNode(to);
  low->upper.push_back(up);
  up->lower.push_back(low);
}

/// TODO: add remove from mapping
void Hasse::RemoveNode(const std::vector<int>& node) {
  RecursiveRemoveNode(node);
}

void Hasse::RemoveArc(const std::vector<int>& from,
                      const std::vector<int>& to) {
  ResetCache();
  auto low = GetNode(from);
  auto up = GetNode(to);
  for (auto it = low->upper.begin(); it != low->upper.end(); it++) {
    if ((*it)->data == to) {
      low->upper.erase(it);
      break;
    }
  }

  for (auto it = up->lower.begin(); it != up->lower.end(); it++) {
    if ((*it)->data == from) {
      up->lower.erase(it);
      break;
    }
  }
}

void Hasse::RecursiveRemoveNode(const std::vector<int>& remove_node) {
  ResetCache();
  std::queue<std::vector<int>> q;
  std::set<std::vector<int>> used;
  q.push(remove_node);
  used.insert(remove_node);
  while (!q.empty()) {
    auto top = q.front();
    q.pop();
    auto node = GetNode(top);
    nodes_with_fixed_rank_[node->rank].erase(node);

    for (const auto& nxt : node->upper) {
      if (!used.count(nxt->data)) {
        used.insert(nxt->data);
        q.push(nxt->data);
      }
    }

    while (!node->lower.empty()) {
      const auto& prev = node->lower.back();
      RemoveArc(prev->data, top);
    }

    node->upper.clear();
    node->lower.clear();
  }

  for (const auto& node : used) {
    mapping_.erase(node);
  }
}

void Hasse::Threshold(std::function<bool(std::vector<int>)> is_good) {
  ResetCache();

  std::vector<int> ranks;
  for (const auto& [rank, _] : nodes_with_fixed_rank_) {
    ranks.push_back(rank);
  }
  reverse(ranks.begin(), ranks.end());
  for (auto rank : ranks) {
    auto nodes = GetNodesWithFixedRank(rank);
    for (auto node : nodes) {
      if (!is_good(node->data)) {
        RemoveNode(node->data);
      }
    }
  }
}

void Hasse::ResetCache() {
  cache_degree_.clear();
  cache_incidence_.clear();
}

bool Hasse::In(const std::vector<int>& a, const std::vector<int>& b) {
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

void Hasse::RecursiveAddNode(const std::vector<int>& add_node) {
  if (mapping_.count(add_node)) {
    return;
  }

  ResetCache();
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
  for (auto& [id, node] : mapping_) {
    if (node->upper.empty()) {
      result.push_back(node->data);
    }
  }
  return result;
}

std::vector<std::vector<int>> Hasse::GetAllElements() {
  std::vector<std::vector<int>> result;
  for (auto& [rank, nodes] : nodes_with_fixed_rank_) {
    for (auto node : nodes) {
      result.push_back(node->data);
    }
  }
  return result;
}

int Hasse::Size() {
  return mapping_.size();
}

std::vector<std::vector<int>>& Hasse::IncidenceMatrix(int p, int k) {
  if (k < p) {
    throw std::runtime_error("k < size of node");
  }
  // p <= k
  if (cache_incidence_.count(std::make_pair(k, p))) {
    return cache_incidence_[std::make_pair(k, p)];
  }

  std::vector<Node*> lower = GetNodesWithFixedRank(p);

  std::vector<Node*> upper = GetNodesWithFixedRank(k);
  std::unordered_map<Node*, int> positions;
  for (size_t i = 0; i < upper.size(); i++) {
    auto up = upper[i];
    positions[up] = i;
  }
  std::vector<std::vector<int>> result(lower.size(),
                                       std::vector<int>(upper.size()));
  for (size_t i = 0; i < lower.size(); i++) {
    auto neighbors = lower[i]->GetAllUpper(k);
    for (auto up : neighbors) {
      auto pos = positions[up];
      assert(pos < upper.size() && upper[pos] == up);
      result[i][pos] = 1;
    }
  }
  cache_incidence_[std::make_pair(k, p)] = std::move(result);
  return cache_incidence_[std::make_pair(k, p)];
}

std::vector<std::vector<double>>& Hasse::DegreeMatrix(int p, int k,
                                                      bool weighted) {
  if (k < p) {
    throw std::runtime_error("k < p");
  }
  if (cache_degree_.count(std::make_pair(k, p))) {
    return cache_degree_[std::make_pair(k, p)];
  }

  std::vector<Node*> upper = GetNodesWithFixedRank(k);

  const std::vector<std::vector<int>>& incidence = IncidenceMatrix(p, k);
  if (incidence.empty()) {
    cache_degree_[std::make_pair(k, p)] = {};
    return cache_degree_[std::make_pair(k, p)];
  }

  int size = incidence.size();
  int inter_size = incidence[0].size();
  std::vector<std::vector<double>> result(size, std::vector<double>(size));
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      if (i == j) {
        continue;
      }
      std::vector<double> weights;
      for (size_t k = 0; k < inter_size; k++) {
        if (incidence[i][k] && incidence[j][k]) {
          if (!weighted) {
            result[i][j] = 1.0;
            break;
          } else {
            weights.push_back(upper[k]->weight);
          }
        }
      }
      if (weighted && !weights.empty()) {
        result[i][j] = *std::min_element(weights.begin(), weights.end());
      }
    }
  }
  cache_degree_[std::make_pair(k, p)] = std::move(result);
  return cache_degree_[std::make_pair(k, p)];
}

int Hasse::Dimension() {
  while (!nodes_with_fixed_rank_.empty()) {
    auto it = std::prev(nodes_with_fixed_rank_.end());
    if (it->second.empty()) {
      nodes_with_fixed_rank_.erase(it);
    } else {
      return it->first;
    }
  }
  assert(nodes_with_fixed_rank_.empty());
  return 0;
}

std::vector<std::pair<int, int>> Hasse::FVector() {
  std::vector<std::pair<int, int>> result;
  for (const auto& [rank, nodes] : nodes_with_fixed_rank_) {
    int size = nodes.size();
    if (size != 0 && rank != -1) {
      result.emplace_back(rank, size);
    }
  }
  return result;
}

int Hasse::TotalCount() {
  int result = 0;
  for (const auto& [rank, nodes] : nodes_with_fixed_rank_) {
    int size = nodes.size();
    if (size != 0 && rank != -1) {
      result += size;
    }
  }
  return result;
}

int Hasse::EulerCharacteristic() {
  int result = 0;
  for (const auto& [rank, nodes] : nodes_with_fixed_rank_) {
    int size = nodes.size();
    if (size != 0 && rank != -1) {
      if (rank % 2 == 0) {
        result += size;
      } else {
        result -= size;
      }
    }
  }
  return result;
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

int Hasse::GetPositionInFixedRank(std::vector<int> node) {
  std::vector<Node*> nodes = GetNodesWithFixedRank(GetNode(node)->rank);
  for (size_t i = 0; i < nodes.size(); i++) {
    if (nodes[i]->data == node) {
      return i;
    }
  }
  throw std::runtime_error("unknown node");
}

int Hasse::BettiNumber(int k) {
  /// TODO: add runtime checks
  /// TODO: check code
  int kernel = 0;
  if (k == 0) {
    kernel = this->nodes_with_fixed_rank_[0].size();
  } else {
    const std::vector<std::vector<int>>& inc = IncidenceMatrix(k - 1, k);
    if (!inc.empty()) {
      kernel = (int)(inc[0].size()) - SNF(inc);
    }
  }
  int image = 0;
  if (nodes_with_fixed_rank_[k + 1].size() != 0) {
    const std::vector<std::vector<int>>& inc = IncidenceMatrix(k, k + 1);
    image = SNF(inc);
  }
  return kernel - image;
}

int Hasse::CalculateSign(const std::vector<int>& subset,
                         const std::vector<int>& set) {
  assert(Hasse::In(subset, set));
  int inv = 0;
  for (size_t i = 0, j = 0; i < subset.size(); i++) {
    while (set[j] != subset[i]) {
      j += 1;
    }
    inv += j - i;
  }
  return inv % 2 == 1 ? -1 : +1;
}

MyMatrixInt Hasse::BoundaryMatrix(int k, int p) {
  if (k <= p) {
    throw std::runtime_error("k <= p");
  }

  std::vector<Node*> lower = GetNodesWithFixedRank(p);

  std::vector<Node*> upper = GetNodesWithFixedRank(k);
  std::unordered_map<Node*, int> positions;
  for (size_t i = 0; i < upper.size(); i++) {
    auto up = upper[i];
    positions[up] = i;
  }
  MyMatrixInt result = MyMatrixInt::Zero(lower.size(), upper.size());
  for (size_t i = 0; i < lower.size(); i++) {
    auto neighbors = lower[i]->GetAllUpper(k);
    for (auto up : neighbors) {
      auto pos = positions[up];
      assert(pos < upper.size() && upper[pos] == up);
      result(i, pos) = CalculateSign(lower[i]->data, upper[pos]->data);
    }
  }
  return result;
}

MyMatrixDouble Hasse::LaplacianMatrix(int k, int p, int q, bool weighted) {
  if (p >= k || k >= q) {
    throw std::runtime_error("Should be p < k < q!");
  }
  assert(p < k && k < q);
  if (!weighted) {
    auto b1 = BoundaryMatrix(k, p);
    auto b2 = BoundaryMatrix(q, k);
    auto result = b1.transpose() * b1 + b2 * b2.transpose();
    return result.cast<double>();
  } else {
    // important lines, can't cast from temporary object
    auto raw_b1 = BoundaryMatrix(k, p);
    auto raw_b2 = BoundaryMatrix(q, k);

    auto b1 = raw_b1.cast<double>();
    auto b2 = raw_b2.cast<double>();
    auto wp = WeightedMatrix(p);
    auto wk = WeightedMatrix(k);
    auto wq = WeightedMatrix(q);
    return b1.transpose() * wp.inverse() * b1 * wk +
           wk.inverse() * b2 * wq * b2.transpose();
  }
}

void Merge(Hasse& current, Hasse& other) {
  for (auto& [key, value] : other.mapping_) {
    current.mapping_[key] = std::move(value);
  }
}

std::vector<std::vector<int>> Hasse::Incidence(const std::vector<int>& node,
                                               int k) {
  int p = GetNode(node)->rank;
  if (k < p) {
    throw std::runtime_error("k < size of node");
  }
  const std::vector<std::vector<int>>& mat = IncidenceMatrix(p, k);
  std::vector<std::vector<int>> result;
  std::vector<Node*> upper = GetNodesWithFixedRank(k);
  std::vector<Node*> lower = GetNodesWithFixedRank(p);

  size_t row = GetPositionInFixedRank(node);

  /// TODO: make bitset or just links
  for (size_t col = 0; col < mat[row].size(); col++) {
    if (mat[row][col]) {
      result.push_back(upper[col]->data);
    }
  }
  return result;
}

int Hasse::IncidenceDegree(const std::vector<int>& node, int k) {
  return Incidence(node, k).size();
}

std::vector<std::vector<int>> Hasse::Adjacency(const std::vector<int>& node,
                                               int k) {
  int p = GetNode(node)->rank;
  std::vector<std::vector<int>> result;
  const std::vector<std::vector<double>>& mat = DegreeMatrix(p, k);

  std::vector<Node*> nodes = GetNodesWithFixedRank(p);

  size_t row = GetPositionInFixedRank(node);

  for (size_t col = 0; col < mat[row].size(); col++) {
    if (mat[row][col]) {
      result.push_back(nodes[col]->data);
    }
  }
  return result;
}

double Hasse::Degree(const std::vector<int>& node, int k, bool weighted) {
  if (!weighted) {
    return Adjacency(node, k).size();
  } else {
    auto adj = Adjacency(node, k);
    if (adj.empty()) {
      return 0.0;
    }
    double sum = 0;
    for (auto v : adj) {
      sum += GetNode(v)->weight;
    }
    return sum;
    // return sum / adj.size();
  }
}

std::vector<double> Hasse::DegreeAll(int p, int k, bool weighted) {
  int n = nodes_with_fixed_rank_[p].size();
  std::vector<double> result;
  std::vector<Node*> nodes = GetNodesWithFixedRank(p);
  for (auto node : nodes) {
    result.emplace_back(Degree(node->data, k, weighted));
  }
  return result;
}

double Hasse::Closeness(std::vector<int> node, int max_rank, bool weighted) {
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

  const std::vector<std::vector<double>>& g =
      DegreeMatrix(k, max_rank, weighted);
  int n = g.size();
  size_t start = GetPositionInFixedRank(node);

  double sum_distances = 0.0;
  int visited = 0;
  std::vector<double> dist(n, n + 1);
  dist[start] = 0;

  if (!weighted) {
    std::queue<int> q;
    q.push(start);
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
  } else {
    std::priority_queue<std::pair<double, int>> q;
    q.push({-dist[start], start});
    while (!q.empty()) {
      auto v = q.top().second;
      auto cur_d = -q.top().first;
      q.pop();
      if (cur_d > dist[v]) {
        continue;
      }
      for (int u = 0; u < n; u++) {
        if (g[v][u] > 0 && dist[u] > dist[v] + g[v][u]) {
          dist[u] = dist[v] + g[v][u];
          q.push({-dist[u], u});
        }
      }
    }
  }

  for (int v = 0; v < n; v++) {
    if (dist[v] < n) {
      sum_distances += dist[v];
      visited += 1;
    }
  }

  /// TODO: now norming by number of visisted nodes
  /// including self
  if (visited == 1) {  // node isolated
    return 0.0;
  }

  double avg_sum_distance = 1.0 * sum_distances / (visited - 1);
  avg_sum_distance *= 1.0 * (n - 1) / (visited - 1);
  return 1 / avg_sum_distance;
}

double Hasse::Betweenness(std::vector<int> node, int max_rank, bool weighted) {
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

  const std::vector<std::vector<double>>& g =
      DegreeMatrix(k, max_rank, weighted);
  int n = g.size();

  const std::vector<Node*>& nodes = GetNodesWithFixedRank(k);

  size_t index_chosen_node = GetPositionInFixedRank(node);

  double sum_distances = 0;

  for (size_t i = 0; i < nodes.size(); i++) {
    auto s = nodes[i];
    if (i == index_chosen_node) {
      continue;
    }
    std::vector<double> dist(n, n);
    std::vector<int> shortest_path_cnt(n, 0);
    std::vector<std::pair<int, int>> shortest_path_through_node_cnt(
        n, {0, 0});  // (cnt, visited)
    dist[i] = 0;
    shortest_path_cnt[i] = 1;
    shortest_path_through_node_cnt[i] = {1, 0};
    if (!weighted) {
      std::queue<int> q;
      q.push(i);
      while (!q.empty()) {
        int v = q.front();
        q.pop();
        if (v == index_chosen_node) {
          shortest_path_through_node_cnt[v].second = true;
        }
        for (int u = 0; u < n; u++) {
          if (!g[v][u]) {
            continue;
          }
          if (dist[u] > dist[v] + 1) {
            dist[u] = dist[v] + 1;
            shortest_path_cnt[u] = shortest_path_cnt[v];
            shortest_path_through_node_cnt[u] =
                shortest_path_through_node_cnt[v];
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
    } else {
      std::priority_queue<std::pair<double, int>> q;
      dist[i] = 0;
      q.push({-dist[i], i});

      while (!q.empty()) {
        auto v = q.top().second;
        auto cur_d = -q.top().first;
        q.pop();
        if (v == index_chosen_node) {
          shortest_path_through_node_cnt[v].second = true;
        }
        if (cur_d > dist[v]) {
          continue;
        }

        for (int u = 0; u < n; u++) {
          if (g[v][u] == 0) {
            continue;
          }
          if (dist[u] > dist[v] + g[v][u]) {
            dist[u] = dist[v] + g[v][u];
            shortest_path_cnt[u] = shortest_path_cnt[v];
            shortest_path_through_node_cnt[u] =
                shortest_path_through_node_cnt[v];
            q.push({-dist[u], u});
          } else if (dist[u] == dist[v] + g[v][u]) {
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
    }
    for (int u = i + 1; u < n; u++) {
      if (u == index_chosen_node) {
        continue;
      }
      int num = shortest_path_through_node_cnt[u].first *
                shortest_path_through_node_cnt[u].second;
      int den = shortest_path_cnt[u];
      if (den == 0) {
        continue;
      }
      sum_distances += 1.0 * num / den;
    }
  }
  /// TODO: norm by n, maybe should norm by visited nodes
  int norm = (n - 1) * (n - 2) / 2;
  return sum_distances / norm;
}

std::vector<std::pair<std::vector<int>, double>> Hasse::ClosenessAll(
    int p, int max_rank, bool weighted) {
  int n = nodes_with_fixed_rank_[p].size();
  std::vector<std::pair<std::vector<int>, double>> result;
  std::vector<Node*> nodes = GetNodesWithFixedRank(p);
  for (auto node : nodes) {
    result.emplace_back(node->data, Closeness(node->data, max_rank, weighted));
  }
  return result;
}

std::vector<std::pair<std::vector<int>, double>> Hasse::BetweennessAll(
    int p, int max_rank, bool weighted) {
  size_t n = nodes_with_fixed_rank_[p].size();
  std::vector<std::pair<std::vector<int>, double>> result;
  std::vector<Node*> nodes = GetNodesWithFixedRank(p);
  int total_threads = std::thread::hardware_concurrency();

  DegreeMatrix(p, max_rank, weighted);

  auto Go = [&](int l, int r,
                std::vector<std::pair<std::vector<int>, double>>& res) {
    if (l >= r) {
      return;
    }
    for (size_t i = l; i < r; i++) {
      res.emplace_back(nodes[i]->data,
                       Betweenness(nodes[i]->data, max_rank, weighted));
    }
  };

  int chunk = std::max((int)n / total_threads, 1);

  std::vector<std::vector<std::pair<std::vector<int>, double>>> part(
      total_threads);
  std::vector<std::thread> threads;

  for (size_t i = 0; i < total_threads; i++) {
    size_t l = i * chunk;
    size_t r = std::min((size_t)(i + 1) * chunk, n);
    std::thread thread(Go, l, r, std::ref(part[i]));
    threads.emplace_back(std::move(thread));
  }

  for (auto& thread : threads) {
    thread.join();
  }

  for (const auto& p : part) {
    for (auto x : p) {
      result.push_back(std::move(x));
    }
  }

  return result;
}

void Hasse::AddFunction(std::string name,
                        std::function<double(std::vector<int>)> func) {
  custom_function_[name] = func;
}

void Hasse::RemoveFunction(std::string name) {
  custom_function_.erase(name);
}

std::vector<std::vector<double>> Hasse::FeaturesMatrix(int rank) {
  auto nodes = GetNodesWithFixedRank(rank);

  std::vector<std::function<double(std::vector<int>)>> features;
  for (const auto& [name, func] : custom_function_) {
    features.emplace_back(func);
  }

  int l = nodes.size();
  int d = custom_function_.size();

  std::vector<std::vector<double>> result(l, std::vector<double>(d));
  for (size_t i = 0; i < l; i++) {
    for (size_t j = 0; j < d; j++) {
      result[i][j] = features[j](nodes[i]->data);
    }
  }
  return result;
}

std::vector<Node*> Hasse::GetNodesWithFixedRank(int rank) {
  std::vector<Node*> nodes(nodes_with_fixed_rank_[rank].begin(),
                           nodes_with_fixed_rank_[rank].end());
  return nodes;
}

void Hasse::ThresholdAbove(std::string name, double threshold) {
  if (!custom_function_.count(name)) {
    throw std::runtime_error("no such function");
  }
  auto& func = custom_function_[name];
  auto is_good = [&](std::vector<int> node) {
    return func(node) < threshold;
  };
  Threshold(is_good);
}

void Hasse::ThresholdBelow(std::string name, double threshold) {
  if (!custom_function_.count(name)) {
    throw std::runtime_error("no such function");
  }
  auto& func = custom_function_[name];
  auto is_good = [&](std::vector<int> node) {
    return func(node) > threshold;
  };
  Threshold(is_good);
}

void Hasse::UpdateWeight(std::vector<int> node, double new_weight) {
  ResetCache();
  assert(is_sorted(node.begin(), node.end()));
  if (!mapping_.count(node)) {
    throw std::runtime_error("no such node exists");
  }
  auto n = mapping_[node].get();
  n->UpdateWeight(new_weight);
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> Hasse::WeightedMatrix(int rank) {
  auto nodes = GetNodesWithFixedRank(rank);
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> result(nodes.size());
  for (size_t i = 0; i < nodes.size(); i++) {
    result.diagonal()[i] = nodes[i]->weight;
  }
  return std::move(result);
}
