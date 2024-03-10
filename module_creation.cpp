#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "src/HyperGraph.h"
#include "src/SimplicialComplex.h"

namespace py = pybind11;

PYBIND11_MODULE(simpl, m) {
    m.doc() = "simplex and etc";

    py::class_<SimplicialComplex>(m, "SimplicialComplex")
        .def(py::init<>())
        .def("AddComplex", &SimplicialComplex::AddComplex)
        .def("RemoveComplex", &SimplicialComplex::RemoveComplex);

    py::class_<HyperGraph>(m, "HyperGraph")
        .def(py::init<>())
        .def("AddEdge", &HyperGraph::AddEdge)
        .def("RemoveEdge", &HyperGraph::RemoveEdge);
}
