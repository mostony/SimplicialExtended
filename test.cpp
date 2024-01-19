#include <pybind11/pybind11.h>

#include <cassert>
#include <iostream>

#include "Hasse.h"
#include "HyperGraph.h"
#include "SimplicialComplex.h"

namespace py = pybind11;

void TestHasse() {
    std::cerr << "TestHasse Started...\n";
    Hasse h;
    h.AddArc({1, 2}, {3, 4, 5});
    h.DebugPrintAll();
    std::cerr << "TestHasse Finished...\n";
}

void TestHyper() {
    std::cerr << "TestHyper Started...\n";

    HyperGraph hyper;
    hyper.AddEdge({1, 2});
    hyper.AddEdge({2, 3});
    hyper.AddEdge({1, 2, 3});
    hyper.RemoveEdge({1, 2});

    hyper.Debug();
    std::cerr << "TestHyper Finished...\n";
}

void TestSimplex() {
    std::cerr << "TestSimplex Started...\n";
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3});
    simpl.RemoveComplex({1, 2});
    // 1 2 3
    // 13 23
    simpl.Debug();
    std::cerr << "TestSimplex Finished...\n";
}

PYBIND11_MODULE(simpl, m) {
    m.doc() = "simplex and etc";
    m.def("TestHasse", &TestHasse, "Testing hasse impl.");
    m.def("TestHyper", &TestHyper, "Testing hypergraph impl.");
    m.def("TestSimplex", &TestSimplex, "Testing simplex impl.");

    // py::class_<Pet>(m, "Pet")
    //     .def(py::init<const std::string &>())  // our constructor
    //     .def("setName", &Pet::setName)         // Expose member methods
    //     .def("getName", &Pet::getName)         // Think about the syntax "&Pet then "::" and the method name
    //     .def_readwrite("name", &Pet::name);    // Expose member variables
}

// int main() {
//     std::cerr << "Start testing...\n";
//     TestHasse();
//     TestHyper();
//     TestSimplex();
// }
