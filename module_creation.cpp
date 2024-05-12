#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"
#include "pybind11/complex.h"
#include <pybind11/eigen.h>

#include "src/CombinatorialComplex.h"
#include "src/Graph.h"
#include "src/HyperGraph.h"
#include "src/SimplicialComplex.h"

namespace py = pybind11;

using namespace pybind11::literals;

PYBIND11_MODULE(simpl, m) {
  m.doc() =
      "Library for SimplicialComplex, Graph, HyperGraph and "
      "CombinatorialComplex";

  py::class_<SimplicialComplex>(m, "SimplicialComplex")
      .def(py::init<>())
      .def("AddFunction", &SimplicialComplex::AddFunction, "name"_a, "func"_a)
      .def("RemoveFunction", &SimplicialComplex::RemoveFunction, "name"_a)
      .def("FeaturesMatrix", &SimplicialComplex::FeaturesMatrix, "rank"_a)
      .def("UpdateWeight", &SimplicialComplex::UpdateWeight, "node"_a, "rank"_a)
      .def("Weights", &SimplicialComplex::Weights, "rank"_a)
      .def("LaplacianMatrix", &SimplicialComplex::LaplacianMatrix, "k"_a, "p"_a,
           "q"_a, "weighted"_a = false)
      .def("EigenValues", &SimplicialComplex::EigenValues, "k"_a, "p"_a, "q"_a,
           "weighted"_a, "cnt"_a)
      .def("EigenValuesAll", &SimplicialComplex::EigenValuesAll, "k"_a, "p"_a,
           "q"_a, "weighted"_a)
      .def("ThresholdAbove", &SimplicialComplex::ThresholdAbove, "name"_a,
           "threshold"_a)
      .def("ThresholdBelow", &SimplicialComplex::ThresholdBelow, "name"_a,
           "threshold"_a)
      .def("AddSimplex", &SimplicialComplex::AddSimplex, "simplex"_a)
      .def("BoundaryMatrix", &SimplicialComplex::BoundaryMatrix, "k"_a, "p"_a)
      .def("RemoveSimplex", &SimplicialComplex::RemoveSimplex, "simplex"_a)
      .def("Incidence", &SimplicialComplex::Incidence, "node"_a, "k"_a)
      .def("IncidenceDegree", &SimplicialComplex::IncidenceDegree, "node"_a,
           "k"_a)
      .def("Adjacency", &SimplicialComplex::Adjacency, "node"_a, "k"_a)
      .def("Degree", &SimplicialComplex::Degree, "node"_a, "k"_a,
           "weighted"_a = false)
      .def("DegreeAll", &SimplicialComplex::DegreeAll, "p"_a, "k"_a,
           "weighted"_a = false)
      .def("BettiNumber", &SimplicialComplex::BettiNumber, "k"_a)
      .def("Closeness", &SimplicialComplex::Closeness, "node"_a, "p"_a,
           "weighted"_a = false)
      .def("ClosenessAll", &SimplicialComplex::ClosenessAll, "k"_a, "p"_a,
           "weighted"_a = false)
      .def("Betweenness", &SimplicialComplex::Betweenness, "node"_a, "p"_a,
           "weighted"_a = false)
      .def("BetweennessAll", &SimplicialComplex::BetweennessAll, "k"_a, "p"_a,
           "weighted"_a = false)
      .def("GetMaxSimplices", &SimplicialComplex::GetMaxSimplices)
      .def("GetElementsWithRank", &SimplicialComplex::GetElementsWithRank)
      .def("GetAll", &SimplicialComplex::GetAll)
      .def("TotalCount", &SimplicialComplex::TotalCount)
      .def("FVector", &SimplicialComplex::FVector)
      .def("Dimension", &SimplicialComplex::Dimension)
      .def("EulerCharacteristic", &SimplicialComplex::EulerCharacteristic)
      .def("BuildFromDowkerComplex", &SimplicialComplex::BuildFromDowkerComplex,
           "binary"_a, "on_column"_a = false)
      .def_static("CreateCliqueGraph", &SimplicialComplex::CreateCliqueGraph,
                  "g"_a, "k"_a, "method"_a, "threads"_a)
      .def("Clear", &SimplicialComplex::Clear);

  py::class_<HyperGraph>(m, "HyperGraph")
      .def(py::init<>())
      .def("AddFunction", &HyperGraph::AddFunction)
      .def("RemoveFunction", &HyperGraph::RemoveFunction)
      .def("FeaturesMatrix", &HyperGraph::FeaturesMatrix)
      .def("UpdateWeight", &HyperGraph::UpdateWeight)
      .def("Weights", &HyperGraph::Weights)
      .def("LaplacianMatrix", &HyperGraph::LaplacianMatrix)
      .def("EigenValues", &HyperGraph::EigenValues)
      .def("EigenValuesAll", &HyperGraph::EigenValuesAll)
      .def("ThresholdAbove", &HyperGraph::ThresholdAbove)
      .def("ThresholdBelow", &HyperGraph::ThresholdBelow)
      .def("AddEdge", &HyperGraph::AddEdge)
      .def("RemoveEdge", &HyperGraph::RemoveEdge)
      .def("BoundaryMatrix", &HyperGraph::BoundaryMatrix)
      .def("Incidence", &HyperGraph::Incidence)
      .def("IncidenceDegree", &HyperGraph::IncidenceDegree)
      .def("Adjacency", &HyperGraph::Adjacency)
      .def("Degree", &HyperGraph::Degree)
      .def("DegreeAll", &HyperGraph::DegreeAll)
      .def("BettiNumber", &HyperGraph::BettiNumber)
      .def("Closeness", &HyperGraph::Closeness)
      .def("ClosenessAll", &HyperGraph::ClosenessAll)
      .def("Betweenness", &HyperGraph::Betweenness)
      .def("BetweennessAll", &HyperGraph::BetweennessAll)
      .def("GetEdges", &HyperGraph::GetEdges)
      .def("GetElementsWithRank", &HyperGraph::GetElementsWithRank)
      .def("GetAll", &HyperGraph::GetAll)
      .def("TotalCount", &HyperGraph::TotalCount)
      .def("FVector", &HyperGraph::FVector)
      .def("Dimension", &HyperGraph::Dimension)
      .def("EulerCharacteristic", &HyperGraph::EulerCharacteristic)
      .def("Clear", &HyperGraph::Clear);

  py::class_<Graph>(m, "Graph")
      .def(py::init<>())
      .def("AddFunction", &Graph::AddFunction)
      .def("RemoveFunction", &Graph::RemoveFunction)
      .def("FeaturesMatrix", &Graph::FeaturesMatrix)
      .def("UpdateWeight", &Graph::UpdateWeight)
      .def("Weights", &Graph::Weights)
      .def("LaplacianMatrix", &Graph::LaplacianMatrix)
      .def("EigenValues", &Graph::EigenValues)
      .def("EigenValuesAll", &Graph::EigenValuesAll)
      .def("ThresholdAbove", &Graph::ThresholdAbove)
      .def("ThresholdBelow", &Graph::ThresholdBelow)
      .def("AddEdge", &Graph::AddEdge)
      .def("RemoveEdge", &Graph::RemoveEdge)
      .def("BoundaryMatrix", &Graph::BoundaryMatrix)
      .def("Incidence", &Graph::Incidence)
      .def("IncidenceDegree", &Graph::IncidenceDegree)
      .def("Adjacency", &Graph::Adjacency)
      .def("Degree", &Graph::Degree)
      .def("DegreeAll", &Graph::DegreeAll)
      .def("BettiNumber", &Graph::BettiNumber)
      .def("Closeness", &Graph::Closeness)
      .def("ClosenessAll", &Graph::ClosenessAll)
      .def("Betweenness", &Graph::Betweenness)
      .def("BetweennessAll", &Graph::BetweennessAll)
      .def("GetEdges", &Graph::GetEdges)
      .def("GetElementsWithRank", &Graph::GetElementsWithRank)
      .def("GetAll", &Graph::GetAll)
      .def("TotalCount", &Graph::TotalCount)
      .def("FVector", &Graph::FVector)
      .def("Dimension", &Graph::Dimension)
      .def("EulerCharacteristic", &Graph::EulerCharacteristic)
      .def("Clear", &Graph::Clear);

  py::class_<CombinatorialComplex>(m, "CombinatorialComplex")
      .def(py::init<>())
      .def("AddFunction", &CombinatorialComplex::AddFunction)
      .def("RemoveFunction", &CombinatorialComplex::RemoveFunction)
      .def("FeaturesMatrix", &CombinatorialComplex::FeaturesMatrix)
      .def("UpdateWeight", &CombinatorialComplex::UpdateWeight)
      .def("Weights", &CombinatorialComplex::Weights)
      .def("LaplacianMatrix", &CombinatorialComplex::LaplacianMatrix)
      .def("EigenValues", &CombinatorialComplex::EigenValues)
      .def("EigenValuesAll", &CombinatorialComplex::EigenValuesAll)
      .def("ThresholdAbove", &CombinatorialComplex::ThresholdAbove)
      .def("ThresholdBelow", &CombinatorialComplex::ThresholdBelow)
      .def("Build", &CombinatorialComplex::Build)
      .def("BoundaryMatrix", &CombinatorialComplex::BoundaryMatrix)
      .def("Incidence", &CombinatorialComplex::Incidence)
      .def("IncidenceDegree", &CombinatorialComplex::IncidenceDegree)
      .def("Adjacency", &CombinatorialComplex::Adjacency)
      .def("Degree", &CombinatorialComplex::Degree)
      .def("DegreeAll", &CombinatorialComplex::DegreeAll)
      .def("BettiNumber", &CombinatorialComplex::BettiNumber)
      .def("Closeness", &CombinatorialComplex::Closeness)
      .def("ClosenessAll", &CombinatorialComplex::ClosenessAll)
      .def("Betweenness", &CombinatorialComplex::Betweenness)
      .def("BetweennessAll", &CombinatorialComplex::BetweennessAll)
      .def("GetEdges", &CombinatorialComplex::GetSubsets)
      .def("GetElementsWithRank", &CombinatorialComplex::GetElementsWithRank)
      .def("GetAll", &CombinatorialComplex::GetAll)
      .def("TotalCount", &CombinatorialComplex::TotalCount)
      .def("FVector", &CombinatorialComplex::FVector)
      .def("Dimension", &CombinatorialComplex::Dimension)
      .def("EulerCharacteristic", &CombinatorialComplex::EulerCharacteristic)
      .def("Clear", &CombinatorialComplex::Clear);
}
