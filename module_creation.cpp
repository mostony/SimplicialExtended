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
      .def("IncidenceDegree", &SimplicialComplex::IncidenceDegree)
      .def("Adjacency", &SimplicialComplex::Adjacency)
      .def("Degree", &SimplicialComplex::Degree)
      .def("BettiNumber", &SimplicialComplex::BettiNumber)
      .def("Closeness", &SimplicialComplex::Closeness)
      .def("Betweenness", &SimplicialComplex::Betweenness)
      .def("GetMaxSimplices", &SimplicialComplex::GetMaxSimplices)
      .def("TotalCount", &SimplicialComplex::TotalCount)
      .def("FVector", &SimplicialComplex::FVector)
      .def("Dimension", &SimplicialComplex::Dimension)
      .def("EulerCharacteristic", &SimplicialComplex::EulerCharacteristic);

  py::class_<HyperGraph>(m, "HyperGraph")
      .def(py::init<>())
      .def("AddEdge", &HyperGraph::AddEdge)
      .def("RemoveEdge", &HyperGraph::RemoveEdge)
      .def("Incidence", &HyperGraph::Incidence)
      .def("IncidenceDegree", &HyperGraph::IncidenceDegree)
      .def("Adjacency", &HyperGraph::Adjacency)
      .def("Degree", &HyperGraph::Degree)
      .def("BettiNumber", &HyperGraph::BettiNumber)
      .def("Closeness", &HyperGraph::Closeness)
      .def("Betweenness", &HyperGraph::Betweenness)
      .def("GetEdges", &HyperGraph::GetEdges)
      .def("TotalCount", &HyperGraph::TotalCount)
      .def("FVector", &HyperGraph::FVector)
      .def("Dimension", &HyperGraph::Dimension)
      .def("EulerCharacteristic", &HyperGraph::EulerCharacteristic);

  py::class_<Graph>(m, "Graph")
      .def(py::init<>())
      .def("AddEdge", &Graph::AddEdge)
      .def("RemoveEdge", &Graph::RemoveEdge)
      .def("Incidence", &Graph::Incidence)
      .def("IncidenceDegree", &Graph::IncidenceDegree)
      .def("Adjacency", &Graph::Adjacency)
      .def("Degree", &Graph::Degree)
      .def("BettiNumber", &Graph::BettiNumber)
      .def("Closeness", &Graph::Closeness)
      .def("GetEdges", &Graph::GetEdges)
      .def("TotalCount", &Graph::TotalCount)
      .def("FVector", &Graph::FVector)
      .def("Dimension", &Graph::Dimension)
      .def("EulerCharacteristic", &Graph::EulerCharacteristic);

  py::class_<CombinatorialComplex>(m, "CombinatorialComplex")
      .def(py::init<>())
      .def("Build", &CombinatorialComplex::Build)
      .def("Incidence", &CombinatorialComplex::Incidence)
      .def("IncidenceDegree", &CombinatorialComplex::IncidenceDegree)
      .def("Adjacency", &CombinatorialComplex::Adjacency)
      .def("Degree", &CombinatorialComplex::Degree)
      .def("BettiNumber", &CombinatorialComplex::BettiNumber)
      .def("Closeness", &CombinatorialComplex::Closeness)
      .def("Betweenness", &CombinatorialComplex::Betweenness)
      .def("GetEdges", &CombinatorialComplex::GetSubsets)
      .def("TotalCount", &CombinatorialComplex::TotalCount)
      .def("FVector", &CombinatorialComplex::FVector)
      .def("Dimension", &CombinatorialComplex::Dimension)
      .def("EulerCharacteristic", &CombinatorialComplex::EulerCharacteristic);
}
