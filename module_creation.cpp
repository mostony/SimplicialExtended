#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "src/CombinatorialComplex.h"
#include "src/Graph.h"
#include "src/HyperGraph.h"
#include "src/SimplicialComplex.h"
namespace py = pybind11;

PYBIND11_MODULE(simpl, m) {
  m.doc() = "SimplicialComplex, HyperGraph & CombinatorialComplex";

  py::class_<SimplicialComplex>(m, "SimplicialComplex")
      .def(py::init<>())
      .def("AddComplex", &SimplicialComplex::AddComplex)
      .def("RemoveComplex", &SimplicialComplex::RemoveComplex)
      .def("Incidence", &SimplicialComplex::Incidence)
      .def("Degree", &SimplicialComplex::Degree)
      .def("BettiNumber", &SimplicialComplex::BettiNumber)
      .def("Closeness", &SimplicialComplex::Closeness)
      .def("Betweenness", &SimplicialComplex::Betweenness)
      .def("GetMaxSimplices", &SimplicialComplex::GetMaxSimplices);

  py::class_<HyperGraph>(m, "HyperGraph")
      .def(py::init<>())
      .def("AddEdge", &HyperGraph::AddEdge)
      .def("RemoveEdge", &HyperGraph::RemoveEdge)
      .def("Incidence", &HyperGraph::Incidence)
      .def("Degree", &HyperGraph::Degree)
      .def("BettiNumber", &HyperGraph::BettiNumber)
      .def("Closeness", &HyperGraph::Closeness)
      .def("Betweenness", &HyperGraph::Betweenness)
      .def("GetEdges", &HyperGraph::GetEdges);

  py::class_<Graph>(m, "Graph")
      .def(py::init<>())
      .def("AddEdge", &Graph::AddEdge)
      .def("RemoveEdge", &Graph::RemoveEdge)
      .def("Incidence", &Graph::Incidence)
      .def("Degree", &Graph::Degree)
      .def("BettiNumber", &Graph::BettiNumber)
      .def("Closeness", &Graph::Closeness)
      .def("GetEdges", &Graph::GetEdges);

  py::class_<CombinatorialComplex>(m, "CombinatorialComplex")
      .def(py::init<>())
      .def("Build", &CombinatorialComplex::Build)
      .def("Incidence", &CombinatorialComplex::Incidence)
      .def("Degree", &CombinatorialComplex::Degree)
      .def("BettiNumber", &CombinatorialComplex::BettiNumber)
      .def("Closeness", &CombinatorialComplex::Closeness)
      .def("Betweenness", &CombinatorialComplex::Betweenness)
      .def("GetEdges", &CombinatorialComplex::GetSubsets);
}
