//
// Created by Anton Mosin on 21.11.2023.
//

#include "Cell.h"

void Cell::Print() {
    std::cout << "id is: " << id_ << "\n";
    std::cout << "depth is: " << depth_ << "\n";
    std::cout << "sons are: ";
    for (auto son: sons_) {
        std::cout << son->id_ << " ";
    }
    std::cout << "\n";

    std::cout << "parents are: ";
    for (auto par: parents_) {
        std::cout << par->id_ << " ";
    }
    std::cout << "\n";
}

int Cell::GetDepth() {
    return depth_;
}

void Cell::AddEdge(std::shared_ptr<Cell> parent, std::shared_ptr<Cell> son) {
    parent->sons_.push_back(son);
    son->parents_.push_back(parent);
    son->depth_ = parent->depth_ + 1;
}

Cell::Cell(int id) {
    id_ = id;
    depth_ = 0;
}

void Cell::RemoveEdge(std::shared_ptr<Cell> parent, std::shared_ptr<Cell> son) {
    // remove son from parent list
    for (auto it = parent->sons_.begin(); it != parent->sons_.end(); it++) {
        if (*it == son) {
            parent->sons_.erase(it);
            break;
        }
    }

    // remove parent from son list
    for (auto it = son->parents_.begin(); it != son->parents_.end(); it++) {
        if (*it == parent) {
            son->parents_.erase(it);
            break;
        }
    }
}

void Cell::RemoveCell(std::shared_ptr<Cell> cell) {
    for (auto &son : cell->sons_) {
        RemoveEdge(cell, son);
    }
    for (auto &par : cell->parents_) {
        RemoveEdge(par, cell);
    }
}
