#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <iostream>
#include <vector>

#include "src/Hasse.h"
#include "src/HyperGraph.h"
#include "src/SimplicialComplex.h"

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

void TestIncidence() {
    std::cerr << "TestIncidence Started...\n";
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3, 4});
    simpl.AddComplex({1, 5, 6});
    for (int k = 1; k <= 4; k++) {
        std::cerr << k << ": ";
        auto incidence_k = simpl.Incidence({1}, k);
        for (auto v : incidence_k) {
            Node::PrintVec(v);
            std::cerr << ", ";
        }
        std::cerr << "\n";
    }
    std::cerr << "TestIncidence Finished...\n";
}

void TestDegree() {
    std::cerr << "TestDegree Started...\n";
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3, 4});
    simpl.AddComplex({1, 5, 6});
    for (int k = 1; k <= 4; k++) {
        std::cerr << k << ": ";
        auto degree_k = simpl.Degree({1}, k);
        for (auto v : degree_k) {
            Node::PrintVec(v);
            std::cerr << ", ";
        }
        std::cerr << "\n";
    }
    std::cerr << "TestDegree Finished...\n";
}

PYBIND11_MODULE(simpl, m) {
    m.doc() = "simplex and etc";
    m.def("TestHasse", &TestHasse, "Testing hasse impl.");
    m.def("TestHyper", &TestHyper, "Testing hypergraph impl.");
    m.def("TestSimplex", &TestSimplex, "Testing simplex impl.");
    m.def("TestIncidence", &TestIncidence, "Testing incidence impl.");
    m.def("TestDegree", &TestDegree, "Testing degree impl.");

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
