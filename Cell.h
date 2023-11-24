//
// Created by Anton Mosin on 21.11.2023.
//

#ifndef SIMPL_CELL_H
#define SIMPL_CELL_H

#include <vector>
#include <memory>
#include <iostream>

class Cell {
public:
    // i guess inheritance there
    explicit Cell(int id);

    void Print();
    int GetDepth();
    static void AddEdge(std::shared_ptr<Cell> parent, std::shared_ptr<Cell> son);
    static void RemoveEdge(std::shared_ptr<Cell> parent, std::shared_ptr<Cell> son);
    static void RemoveCell(std::shared_ptr<Cell> cell);
private:
    int depth_;
    std::vector<std::shared_ptr<Cell>> sons_, parents_;
    std::vector<int> data_;
    int id_;
    int rank_;
};

#endif //SIMPL_CELL_H
