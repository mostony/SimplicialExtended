#pragma once

#include <algorithm>
#include <exception>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "Hasse.h"

class CombinatorialComplex {
   public:
    int HasseSize();

    void Build(std::vector<std::vector<int>> data, std::vector<int> rank) {}

    std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

    std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

   private:
    Hasse hasse_;
};
