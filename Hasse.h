//
// Created by Anton Mosin on 21.11.2023.
//

#ifndef SIMPL_HASSE_H
#define SIMPL_HASSE_H

#include "Cell.h"
#include <vector>
#include <memory>

class Hasse {
public:
    void PrintAll();
    void AddEdge(std::shared_ptr<Cell> a, std::shared_ptr<Cell> b);
    void AddCell(std::shared_ptr<Cell> cell);
    void RemoveCell(std::shared_ptr<Cell> cell);
private:
    std::vector<std::shared_ptr<Cell>> cells_;
};

#endif //SIMPL_HASSE_H
