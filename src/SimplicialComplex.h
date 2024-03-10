#pragma once

#include <algorithm>
#include <exception>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "Hasse.h"

// TODO: inheritance from Hasse
class SimplicialComplex {
   public:
    void AddComplex(std::vector<int> complex);

    void RemoveComplex(std::vector<int> complex);

    int HasseSize();

    /// merge two non intersecting SimplicialComplex
    friend void Merge(SimplicialComplex* current, SimplicialComplex* other);

    // TODO: change method to enum
    static SimplicialComplex* CreateCliqueGraph(
        const std::vector<std::vector<int>>& g, int k, int method = 0,
        int total_threads = -1);

    std::vector<std::vector<int>> GetMaxSimplices();
    void AddArc(const std::vector<int>& from, const std::vector<int>& to);

    std::vector<std::vector<int>> Incidence(std::vector<int> node, int k);

    // TODO: maybe remove self node from degree???
    std::vector<std::vector<int>> Degree(std::vector<int> node, int k);

   private:
    Hasse hasse_;
};
