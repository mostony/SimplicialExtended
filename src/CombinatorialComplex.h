#pragma once

#include <vector>

#include "Hasse.h"

class CombinatorialComplex {
   public:
    int HasseSize();

    void Build(const std::vector<std::vector<int>>& data);

    std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

    std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

    int BettiNumber(int k);

    double Closeness(std::vector<int> node, int max_rank);

   private:
    Hasse hasse_;
};
