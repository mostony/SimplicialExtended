#pragma once

#include <iostream>
#include <vector>

struct Node {
    int depth = 0;
    // TODO change on Node*
    std::vector<std::vector<int>> upper, lower;
    std::vector<int> data;
    int rank = 0;

    Node() {}

    Node(const std::vector<int>& node) { data = node; }

    static void PrintVec(const std::vector<int>& vec) {
        std::cerr << "{";
        for (int i = 0; i < vec.size(); i++) {
            if (i) {
                std::cerr << ", ";
            }
            std::cerr << vec[i];
        }
        std::cerr << "}";
    }

    void Debug() {
        std::cerr << "Node : ";
        PrintVec(data);
        std::cerr << "\n";

        std::cerr << "sons : ";
        for (auto x : upper) {
            PrintVec(x);
            std::cerr << " ";
        }
        std::cerr << "\n";

        std::cerr << "lower : ";
        for (auto x : lower) {
            PrintVec(x);
            std::cerr << " ";
        }
        std::cerr << "\n";
    }
};
