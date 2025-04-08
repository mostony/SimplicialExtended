#pragma once

#include <vector>

struct Node {
    Node();

    Node(const std::vector<int>& node);

    Node(const std::vector<int>& node, int init_rank);

    void UpdateWeight(double new_weight);

    std::vector<Node*> GetAllUpper(int upper_rank);

    std::vector<Node*> upper, lower;
    std::vector<int> data;
    int rank = 0;
    int size = 0;
    double weight = 1;
};
