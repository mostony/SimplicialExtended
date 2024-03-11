#include "CombinatorialComplex.h"

int CombinatorialComplex::HasseSize() {
    return hasse_.Size();
}

void CombinatorialComplex::Build(const std::vector<std::vector<int>>& data) {
    std::vector<Node*> nodes;
    for (size_t i = 0; i < data.size(); i++) {
        nodes.push_back(hasse_.GetNode(data[i]));
    }

    /// check if a \in b
    auto In = [](const std::vector<int>& a, const std::vector<int>& b) -> bool {
        // TODO: make O(n)
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
    };

    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            if (i == j) {
                continue;
            }
            if (In(data[i], data[j])) {
                hasse_.AddArc(data[i], data[j]);
            }
        }
    }
}

std::vector<std::vector<int>> CombinatorialComplex::Incidence(
    std::vector<int> node, int k) {
    return hasse_.Incidence(node, k);
}

std::vector<std::vector<int>> CombinatorialComplex::Degree(
    std::vector<int> node, int k) {
    return hasse_.Degree(node, k);
}
