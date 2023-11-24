//
// Created by Anton Mosin on 21.11.2023.
//
#include "Hasse.h"
#include <iostream>

int main() {
    // simple test just for test working
    Hasse hasse;

    std::shared_ptr<Cell> c1(new Cell(0));
    std::shared_ptr<Cell> c2(new Cell(1));
    hasse.AddCell(c1);
    hasse.AddCell(c2);
    hasse.AddEdge(c1, c2);

    hasse.PrintAll();
    std::cout << "\n";

    hasse.RemoveCell(c2);

    hasse.PrintAll();
}