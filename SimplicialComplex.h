#pragma once
#include <map>
#include <vector>

#include "Hasse.h"

class SimplicialComplex {
    typedef int32_t VertexId;

   public:
    void AddComplex(std::vector<VertexId> complex);

    void RemoveComplex(std::vector<VertexId> complex);

    void Debug();

   private:
    Hasse hasse_;
};
