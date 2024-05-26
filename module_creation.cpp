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
      "Library for creating complex networks and calculating features on them. "
      "Supported networks: simplicial complexes, "
      "graphs, hypergraphs and "
      "combinatorial complexes";

  py::class_<SimplicialComplex>(m, "SimplicialComplex")
      .def(py::init<>(), "Creating a new simplicial complex")

      .def("AddSimplex", &SimplicialComplex::AddSimplex, "simplex"_a,
           "Add a simplex to the simplicial complex")
      .def("RemoveSimplex", &SimplicialComplex::RemoveSimplex, "simplex"_a,
           "Remove a simplex from the simplicial complex")

      .def("AddFunction", &SimplicialComplex::AddFunction, "name"_a, "func"_a,
           "Add a function on simplices with given name")
      .def("RemoveFunction", &SimplicialComplex::RemoveFunction, "name"_a,
           "Remove a function from simplicial complex by given name")
      .def("FeaturesMatrix", &SimplicialComplex::FeaturesMatrix, "dim"_a,
           "Get the feature matrix with simplices on given dimension")
      .def("UpdateWeight", &SimplicialComplex::UpdateWeight, "node"_a,
           "weight"_a, "Update the weight of a simplex")
      .def("Weights", &SimplicialComplex::Weights, "rank"_a,
           "Get the weights of all nodes at a given rank")
      .def("ThresholdAbove", &SimplicialComplex::ThresholdAbove, "name"_a,
           "threshold"_a,
           "Threshold the simplicial complex by removing elements on which "
           "function with given name higher or equal threshold value. Function "
           "should me monotonically non-decreasing")
      .def("ThresholdBelow", &SimplicialComplex::ThresholdBelow, "name"_a,
           "threshold"_a,
           "Threshold the simplicial complex by removing elements on which "
           "function with given name less or equal threshold value. Function "
           "should me monotonically non-increasing")

      .def("BoundaryMatrix", &SimplicialComplex::BoundaryMatrix, "k"_a, "p"_a,
           "Get the boundary matrix. Should be: k > p")
      .def("Incidence", &SimplicialComplex::Incidence, "node"_a, "k"_a,
           "Get all nodes k-incidence to given node")
      .def("IncidenceDegree", &SimplicialComplex::IncidenceDegree, "node"_a,
           "k"_a, "Get number of nodes k-incidence to given node")
      .def("Adjacency", &SimplicialComplex::Adjacency, "node"_a, "k"_a,
           "Get all nodes k-adjacence to given node")
      .def("Degree", &SimplicialComplex::Degree, "node"_a, "k"_a,
           "weighted"_a = false,
           "Get sum of weights of all k-adjacence to given node.")
      .def("DegreeAll", &SimplicialComplex::DegreeAll, "p"_a, "k"_a,
           "weighted"_a = false,
           "Get vector of degrees of all nodes with rank=p with k-adjacence")

      .def("LaplacianMatrix", &SimplicialComplex::LaplacianMatrix, "k"_a, "p"_a,
           "q"_a, "weighted"_a = false,
           "Get the laplacian matrix. Should be: p < k < q")
      .def("EigenValues", &SimplicialComplex::EigenValues, "k"_a, "p"_a, "q"_a,
           "weighted"_a, "cnt"_a,
           "Get the cnt eigen values of laplacian matrix. Should be: p < k < q")
      .def("EigenValuesAll", &SimplicialComplex::EigenValuesAll, "k"_a, "p"_a,
           "q"_a, "weighted"_a,
           "Get all eigen values of laplacian matrix. Should be: p < k < q")

      .def("BettiNumber", &SimplicialComplex::BettiNumber, "k"_a,
           "Get the k betti number of simplicial complex")
      .def("Closeness", &SimplicialComplex::Closeness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate closeness for given node via nodes with rank=q")
      .def("ClosenessAll", &SimplicialComplex::ClosenessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate closeness for all node with rank=k via nodes with rank=q")
      .def("Betweenness", &SimplicialComplex::Betweenness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate betweenness for given node via nodes with rank=q")
      .def("BetweennessAll", &SimplicialComplex::BetweennessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate betweenness for all node with rank=k via nodes with "
           "rank=q")

      .def("GetMaxSimplices", &SimplicialComplex::GetMaxSimplices,
           "Get all maximum simplices of complex")
      .def("GetSimplicesWithDimension", &SimplicialComplex::GetElementsWithRank,
           "dim"_a, "Get all simplices with given dimension")
      .def("GetAllSimplices", &SimplicialComplex::GetAll, "Get all simplices")
      .def("TotalCountOfSimplices", &SimplicialComplex::TotalCount,
           "Get number of simplices in complex")
      .def("FVector", &SimplicialComplex::FVector,
           "Get the f-vector of the simplicial complex")
      .def("Dimension", &SimplicialComplex::Dimension,
           "Get the dimension of the simplicial complex - maximum dimension of "
           "simplices")
      .def("EulerCharacteristic", &SimplicialComplex::EulerCharacteristic,
           "Get the Euler characteristic of the simplicial complex")

      .def("BuildDowkerComplex", &SimplicialComplex::BuildFromBinary,
           "binary"_a, "on_column"_a = false,
           "Build a Dowker complex from a binary matrix with flag on_column. "
           "If on_column=true then simplices is on column, otherwise on rows")
      .def_static("CreateCliqueGraph", &SimplicialComplex::CreateCliqueGraph,
                  "g"_a, "k"_a, "method"_a = 0, "threads"_a = 2,
                  "Create simplicial complex from given graph using cliques "
                  "with size <= k")

      .def("Clear", &SimplicialComplex::Clear, "Clear the simplicial complex");

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
      .def("BuildWithRank", &CombinatorialComplex::BuildWithRank)
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
