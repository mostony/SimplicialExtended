#pragma once

#include <map>

#include "Hasse.h"

class HyperGraph {
   public:
    typedef int32_t Vertex;

    void AddEdge(const std::vector<Vertex>& edge);

    void RemoveEdge(const std::vector<Vertex>& edge);

    void Debug();

   private:
    Hasse hasse_;
};