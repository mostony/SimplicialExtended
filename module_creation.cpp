#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

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

    py::class_<SimplicialComplex>(m, "SimplicialComplex")
        .def(py::init<>())
        .def("AddComplex", &SimplicialComplex::AddComplex)
        .def("RemoveComplex", &SimplicialComplex::RemoveComplex)
        .def("Debug", &SimplicialComplex::Debug);

    py::class_<HyperGraph>(m, "HyperGraph")
        .def(py::init<>())
        .def("AddEdge", &HyperGraph::AddEdge)
        .def("RemoveEdge", &HyperGraph::RemoveEdge)
        .def("Debug", &HyperGraph::Debug);
}
