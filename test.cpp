#include "Hasse.h"
#include "HyperGraph.h"
#include <iostream>
#include <cassert>

#include "SimplicialComplex.h"

void TestHasse() {
    Hasse h;
    h.AddArc({1, 2}, {3, 4, 5});
    h.DebugPrintAll();
}

void TestHyper() {
    HyperGraph hyper;
    hyper.AddEdge({1, 2});
    hyper.AddEdge({2, 3});
    hyper.AddEdge({1, 2, 3});
    hyper.RemoveEdge({1, 2});

    hyper.Debug();
}

void TestSimplex() {
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3});
    simpl.RemoveComplex({1, 2});
    // 1 2 3
    // 13 23
    simpl.Debug();
}

int main() {
    std::cerr << "Start testing...\n";
//    TestHasse();
//    TestHyper();
    TestSimplex();

}