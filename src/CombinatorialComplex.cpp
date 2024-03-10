#include "CombinatorialComplex.h"

int CombinatorialComplex::HasseSize() {
    return hasse_.Size();
}

std::vector<std::vector<int>> CombinatorialComplex::Incidence(
    std::vector<int> node, int k) {}

std::vector<std::vector<int>> CombinatorialComplex::Degree(
    std::vector<int> node, int k) {}
