#pragma once
#include <exception>
#include <map>
#include <vector>

#include "Hasse.h"

class SimplicialComplex {
    typedef int32_t VertexId;

   public:
    void AddComplex(std::vector<VertexId> complex);

    void RemoveComplex(std::vector<VertexId> complex);

    void Debug();

    // TODO: change method to enum
    static SimplicialComplex CreateCliqueGraph(const std::vector<std::vector<int>> &g,
                                               int k, int method = 0);

    std::vector<std::vector<int>> GetMaxFaces();

   private:
    Hasse hasse_;
};
