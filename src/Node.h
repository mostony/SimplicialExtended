#pragma once

#include <vector>

struct Node {
  std::vector<std::vector<int>> upper, lower;
  std::vector<int> data;
  int rank = 0;
  int size = 0;
  Node();

  Node(const std::vector<int> &node);
};
