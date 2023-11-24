//
// Created by Anton Mosin on 21.11.2023.
//
#include "Hasse.h"

void Hasse::PrintAll() {
    for (const auto &cell : cells_) {
        cell->Print();
    }
}

void Hasse::AddEdge(std::shared_ptr<Cell> a, std::shared_ptr<Cell> b) {
    Cell::AddEdge(a, b);
}

void Hasse::AddCell(std::shared_ptr<Cell> cell) {
    cells_.push_back(cell);
}

void Hasse::RemoveCell(std::shared_ptr<Cell> cell) {
    // stupid, can better
    Cell::RemoveCell(cell);

    for (auto it = cells_.begin(); it != cells_.end(); it++) {
        if (*it == cell) {
            cells_.erase(it);
            break;
        }
    }
}
