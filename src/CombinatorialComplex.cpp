#include "CombinatorialComplex.h"

#include <cassert>
#include <stdexcept>

void CombinatorialComplex::Build(std::vector<std::vector<int>> data) {
    std::vector<Node*> nodes;
    for (size_t i = 0; i < data.size(); i++) {
        std::sort(data[i].begin(), data[i].end());
        nodes.push_back(hasse_.GetNode(data[i]));
    }

    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            if (i == j) {
                continue;
            }
            if (In(data[i], data[j])) {
                assert(nodes[i]->rank <= nodes[j]->rank);
                hasse_.AddArc(data[i], data[j]);
            }
        }
    }
}

void CombinatorialComplex::BuildWithRank(std::vector<std::vector<int>> data,
                                         std::vector<int> ranks) {
    if (data.size() != ranks.size()) {
        throw std::runtime_error("different sizes of data and ranks");
    }
    std::vector<Node*> nodes;
    for (size_t i = 0; i < data.size(); i++) {
        std::sort(data[i].begin(), data[i].end());
        hasse_.CreateNode(data[i], ranks[i]);
        nodes.push_back(hasse_.GetNode(data[i]));
    }

    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            if (i == j) {
                continue;
            }
            if (In(data[i], data[j])) {
                assert(nodes[i]->rank <= nodes[j]->rank);
                hasse_.AddArc(data[i], data[j]);
            }
        }
    }
}

std::vector<std::vector<int>> CombinatorialComplex::GetSubsets() {
    return GetMaxFaces();
}

void CombinatorialComplex::BuildFromBinary(std::vector<std::vector<int>> binary, bool on_column) {
    Clear();
    if (binary.empty()) {
        return;
    }
    std::vector<std::vector<int>> data;

    if (on_column == false) {
        for (size_t i = 0; i < binary.size(); i++) {
            std::vector<int> edge;
            for (size_t j = 0; j < binary[i].size(); j++) {
                if (binary[i][j]) {
                    edge.push_back(j);
                }
            }
            data.push_back(edge);
        }
    } else {
        for (size_t j = 0; j < binary[0].size(); j++) {
            std::vector<int> edge;
            for (size_t i = 0; i < binary.size(); i++) {
                if (binary[i][j]) {
                    edge.push_back(i);
                }
            }
            data.push_back(edge);
        }
    }
    Build(data);
}
