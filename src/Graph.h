#pragma once

#include "Hasse.h"

class Graph {
   public:
    void AddEdge(int v, int u);

    void RemoveEdge(int v, int u);

   private:
    Hasse hasse_;
};