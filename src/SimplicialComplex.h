#pragma once

#include <algorithm>
#include <exception>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#include "Hasse.h"

class SimplicialComplex {
    typedef int32_t VertexId;

   public:
    void AddComplex(std::vector<VertexId> complex);

    void RemoveComplex(std::vector<VertexId> complex);

    void Debug();

    int HasseSize();

    friend void Merge(SimplicialComplex *current, SimplicialComplex *other);

    // TODO: change method to enum
    static SimplicialComplex *CreateCliqueGraph(
        const std::vector<std::vector<int>> &g, int k, int method = 0,
        int total_threads = -1);

    std::vector<std::vector<int>> GetMaxFaces();
    void AddArc(const std::vector<int> &from, const std::vector<int> &to);

   private:
    Hasse hasse_;
};
