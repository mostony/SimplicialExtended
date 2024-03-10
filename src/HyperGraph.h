#pragma once

#include <map>

#include "Hasse.h"

class HyperGraph {
   public:
    void AddEdge(const std::vector<int>& edge);

    void RemoveEdge(const std::vector<int>& edge);

   private:
    Hasse hasse_;
};