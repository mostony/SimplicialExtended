#pragma once

#include <iostream>
#include <memory>
#include <vector>

struct Node {
    int depth = 0;
    std::vector<std::vector<int>> sons,
        parents;
    std::vector<int> data;
    int rank = 0;

    Node() {}

    Node(const std::vector<int> &node) {
        data = node;
    }

    static void PrintVec(const std::vector<int> &vec) {
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
        for (auto x : sons) {
            PrintVec(x);
            std::cerr << " ";
        }
        std::cerr << "\n";

        std::cerr << "parents : ";
        for (auto x : parents) {
            PrintVec(x);
            std::cerr << " ";
        }
        std::cerr << "\n";
    }
};
