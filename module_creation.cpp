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
      .def(py::init<>(), "Creating a new hypergraph")

      .def("AddFunction", &HyperGraph::AddFunction, "name"_a, "func"_a,
           "Add a function to the hypergraph with the given name")
      .def("RemoveFunction", &HyperGraph::RemoveFunction, "name"_a,
           "Remove a function from the hypergraph by given name")
      .def("FeaturesMatrix", &HyperGraph::FeaturesMatrix, "dim"_a,
           "Get the feature matrix of the hypergraph for the specified "
           "dimension")
      .def("UpdateWeight", &HyperGraph::UpdateWeight, "node"_a, "weight"_a,
           "Update the weight of a hyperedge in the hypergraph")
      .def("Weights", &HyperGraph::Weights, "rank"_a,
           "Get the weights of all nodes in the hypergraph with the specified "
           "rank")
      .def("LaplacianMatrix", &HyperGraph::LaplacianMatrix, "k"_a, "p"_a, "q"_a,
           "weighted"_a = false,
           "Get the Laplacian matrix of the hypergraph. Should be: p < k < q")
      .def("EigenValues", &HyperGraph::EigenValues, "k"_a, "p"_a, "q"_a,
           "weighted"_a, "cnt"_a,
           "Get the cnt eigenvalues of the Laplacian matrix. Should be: p < k "
           "< q")
      .def("EigenValuesAll", &HyperGraph::EigenValuesAll, "k"_a, "p"_a, "q"_a,
           "weighted"_a,
           "Get all the eigenvalues of the Laplacian matrix. Should be: p < k "
           "< q")
      .def("ThresholdAbove", &HyperGraph::ThresholdAbove, "name"_a,
           "threshold"_a,
           "Threshold the hypergraph by retaining elements above the given "
           "threshold for the specified function name. The function should be "
           "monotonically non-decreasing")
      .def("ThresholdBelow", &HyperGraph::ThresholdBelow, "name"_a,
           "threshold"_a,
           "Threshold the hypergraph by retaining elements below the given "
           "threshold for the specified function name. The function should be "
           "monotonically non-increasing")
      .def("AddEdge", &HyperGraph::AddEdge, "edge"_a,
           "Add an hyperedge to the hypergraph")
      .def("RemoveEdge", &HyperGraph::RemoveEdge, "edge"_a,
           "Remove an hyperedge from the hypergraph")
      .def("BoundaryMatrix", &HyperGraph::BoundaryMatrix, "k"_a, "p"_a,
           "Get the boundary matrix of the hypergraph. Should be: k > p")
      .def("Incidence", &HyperGraph::Incidence, "node"_a, "k"_a,
           "Get all nodes with k-incidence to the given node")
      .def("IncidenceDegree", &HyperGraph::IncidenceDegree, "node"_a, "k"_a,
           "Get the number of nodes with k-incidence to the given node")
      .def("Adjacency", &HyperGraph::Adjacency, "node"_a, "k"_a,
           "Get all nodes with k-adjacency to the given node")
      .def("Degree", &HyperGraph::Degree, "node"_a, "k"_a, "weighted"_a = false,
           "Get the sum of weights of all k-adjacent nodes to the given node")
      .def("DegreeAll", &HyperGraph::DegreeAll, "rank"_a, "k"_a,
           "weighted"_a = false,
           "Get the degrees of all nodes with the specified rank and "
           "k-adjacency")
      .def("BettiNumber", &HyperGraph::BettiNumber, "k"_a,
           "Get the k-th Betti number of the hypergraph")
      .def("Closeness", &HyperGraph::Closeness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for the given node via nodes "
           "with rank=q")
      .def("ClosenessAll", &HyperGraph::ClosenessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("Betweenness", &HyperGraph::Betweenness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the betweenness centrality for the given node via nodes "
           "with rank=q")
      .def("BetweennessAll", &HyperGraph::BetweennessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the betweenness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("GetEdges", &HyperGraph::GetEdges, "Get all edges in the hypergraph")
      .def("GetElementsWithRank", &HyperGraph::GetElementsWithRank, "rank"_a,
           "Get all elements of the hypergraph with the specified rank")
      .def("GetAll", &HyperGraph::GetAll, "Get all elements in the hypergraph")
      .def("TotalCount", &HyperGraph::TotalCount,
           "Get the total count of elements in the hypergraph")
      .def("FVector", &HyperGraph::FVector,
           "Get the F-vector of the hypergraph")
      .def("Dimension", &HyperGraph::Dimension,
           "Get the dimension of the hypergraph")
      .def("EulerCharacteristic", &HyperGraph::EulerCharacteristic,
           "Get the Euler characteristic of the hypergraph")

      .def("BuildFromBinary", &HyperGraph::BuildFromBinary, "binary"_a,
           "on_column"_a = false,
           "Build from a binary matrix with flag on_column. "
           "If on_column=true then edges is on column, otherwise on rows")

      .def("Clear", &HyperGraph::Clear, "Clear all elements in the hypergraph");

  py::class_<Graph>(m, "Graph")
      .def(py::init<>(), "Creating a new graph")

      .def(
          "AddFunction", &Graph::AddFunction, "name"_a, "func"_a,
          "Add a function to the graph with the given name and function object")
      .def("RemoveFunction", &Graph::RemoveFunction, "name"_a,
           "Remove a function from the graph by given name")
      .def("FeaturesMatrix", &Graph::FeaturesMatrix, "dim"_a,
           "Get the feature matrix of the graph for the specified dimension")
      .def("UpdateWeight", &Graph::UpdateWeight, "node"_a, "weight"_a,
           "Update the weight of a node/edge in the graph")
      .def("Weights", &Graph::Weights, "rank"_a,
           "Get the weights of all node/edges in the graph with the specified "
           "rank")
      .def("LaplacianMatrix", &Graph::LaplacianMatrix, "k"_a, "p"_a, "q"_a,
           "weighted"_a = false, "Get the Laplacian matrix of the graph")
      .def("EigenValues", &Graph::EigenValues, "k"_a, "p"_a, "q"_a,
           "weighted"_a, "cnt"_a,
           "Get the cnt eigenvalues of the Laplacian matrix")
      .def("EigenValuesAll", &Graph::EigenValuesAll, "k"_a, "p"_a, "q"_a,
           "weighted"_a, "Get all the eigenvalues of the Laplacian matrix")
      .def("ThresholdAbove", &Graph::ThresholdAbove, "name"_a, "threshold"_a,
           "Threshold the graph by retaining elements above the given "
           "threshold for the specified function name. The function should be "
           "monotonically non-decreasing")
      .def("ThresholdBelow", &Graph::ThresholdBelow, "name"_a, "threshold"_a,
           "Threshold the graph by retaining elements below the given "
           "threshold for the specified function name. The function should be "
           "monotonically non-increasing")
      .def("AddEdge", &Graph::AddEdge, "v"_a, "u"_a, "Add an edge to the graph")
      .def("RemoveEdge", &Graph::RemoveEdge, "v"_a, "u"_a,
           "Remove an edge from the graph")
      .def("BoundaryMatrix", &Graph::BoundaryMatrix, "k"_a, "p"_a,
           "Get the boundary matrix of the graph. Should be: k > p")
      .def("Incidence", &Graph::Incidence, "node"_a, "k"_a,
           "Get all nodes with k-incidence to the given node")
      .def("IncidenceDegree", &Graph::IncidenceDegree, "node"_a, "k"_a,
           "Get the number of nodes with k-incidence to the given node")
      .def("Adjacency", &Graph::Adjacency, "node"_a, "k"_a,
           "Get all nodes with k-adjacency to the given node")
      .def("Degree", &Graph::Degree, "node"_a, "k"_a, "weighted"_a = false,
           "Get the sum of weights of all k-adjacent nodes to the given node")
      .def("DegreeAll", &Graph::DegreeAll, "rank"_a, "k"_a,
           "weighted"_a = false,
           "Get the degrees of all nodes with the specified rank and "
           "k-adjacency")
      .def("BettiNumber", &Graph::BettiNumber, "k"_a,
           "Get the k-th Betti number of the graph")
      .def("Closeness", &Graph::Closeness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for the given node via nodes "
           "with rank=q")
      .def("ClosenessAll", &Graph::ClosenessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("Betweenness", &Graph::Betweenness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the betweenness centrality for the given node via nodes "
           "with rank=q")
      .def("BetweennessAll", &Graph::BetweennessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the betweenness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("GetEdges", &Graph::GetEdges, "Get all edges in the graph")
      .def("GetElementsWithRank", &Graph::GetElementsWithRank, "rank"_a,
           "Get all elements of the graph with the specified rank")
      .def("GetAll", &Graph::GetAll, "Get all elements in the graph")
      .def("TotalCount", &Graph::TotalCount,
           "Get the total count of elements in the graph")
      .def("FVector", &Graph::FVector, "Get the F-vector of the graph")
      .def("Dimension", &Graph::Dimension, "Get the dimension of the graph")
      .def("EulerCharacteristic", &Graph::EulerCharacteristic,
           "Get the Euler characteristic of the graph")

      .def(
          "BuildFromBinary", &Graph::BuildFromBinary, "binary"_a,
          "on_column"_a = false,
          "Build graph using cliques from a binary matrix with flag on_column. "
          "If on_column=true then edges is on column, otherwise on rows")

      .def("Clear", &Graph::Clear, "Clear all elements in the graph");

  py::class_<CombinatorialComplex>(m, "CombinatorialComplex")
      .def(py::init<>(), "Creating a new combinatorial complex")

      .def("AddFunction", &CombinatorialComplex::AddFunction, "name"_a,
           "func"_a,
           "Add a function to the combinatorial complex with the given name "
           "and function object")
      .def("RemoveFunction", &CombinatorialComplex::RemoveFunction, "name"_a,
           "Remove a function from the combinatorial complex by given name")
      .def("FeaturesMatrix", &CombinatorialComplex::FeaturesMatrix, "dim"_a,
           "Get the feature matrix of the combinatorial complex for the "
           "specified dimension")
      .def("UpdateWeight", &CombinatorialComplex::UpdateWeight, "node"_a,
           "weight"_a,
           "Update the weight of a node in the combinatorial complex")
      .def("Weights", &CombinatorialComplex::Weights, "rank"_a,
           "Get the weights of all nodes in the combinatorial complex with the "
           "specified rank")
      .def("LaplacianMatrix", &CombinatorialComplex::LaplacianMatrix, "k"_a,
           "p"_a, "q"_a, "weighted"_a = false,
           "Get the Laplacian matrix of the combinatorial complex. Should be: "
           "p < k < q")
      .def("EigenValues", &CombinatorialComplex::EigenValues, "k"_a, "p"_a,
           "q"_a, "weighted"_a, "cnt"_a,
           "Get the cnt eigenvalues of the Laplacian matrix. Should be: p < k "
           "< q")
      .def("EigenValuesAll", &CombinatorialComplex::EigenValuesAll, "k"_a,
           "p"_a, "q"_a, "weighted"_a,
           "Get all the eigenvalues of the Laplacian matrix. Should be: p < k "
           "< q")
      .def("ThresholdAbove", &CombinatorialComplex::ThresholdAbove, "name"_a,
           "threshold"_a,
           "Threshold the combinatorial complex by retaining elements above "
           "the given threshold for the specified function name. The function "
           "should be monotonically non-decreasing")
      .def("ThresholdBelow", &CombinatorialComplex::ThresholdBelow, "name"_a,
           "threshold"_a,
           "Threshold the combinatorial complex by retaining elements below "
           "the given threshold for the specified function name. The function "
           "should be monotonically non-increasing")
      .def("Build", &CombinatorialComplex::Build, "elements"_a,
           "Build the combinatorial complex from the given elements")
      .def("BuildWithRank", &CombinatorialComplex::BuildWithRank, "elements"_a,
           "rank"_a,
           "Build the combinatorial complex with the specified rank from the "
           "given elements")
      .def("BoundaryMatrix", &CombinatorialComplex::BoundaryMatrix, "k"_a,
           "p"_a,
           "Get the boundary matrix of the combinatorial complex. Should be: k "
           "> p")
      .def("Incidence", &CombinatorialComplex::Incidence, "node"_a, "k"_a,
           "Get all nodes with k-incidence to the given node")
      .def("IncidenceDegree", &CombinatorialComplex::IncidenceDegree, "node"_a,
           "k"_a, "Get the number of nodes with k-incidence to the given node")
      .def("Adjacency", &CombinatorialComplex::Adjacency, "node"_a, "k"_a,
           "Get all nodes with k-adjacency to the given node")
      .def("Degree", &CombinatorialComplex::Degree, "node"_a, "k"_a,
           "weighted"_a = false,
           "Get the sum of weights of all k-adjacent nodes to the given node")
      .def("DegreeAll", &CombinatorialComplex::DegreeAll, "rank"_a, "k"_a,
           "weighted"_a = false,
           "Get the degrees of all nodes with the specified rank and "
           "k-adjacency")
      .def("BettiNumber", &CombinatorialComplex::BettiNumber, "k"_a,
           "Get the k-th Betti number of the combinatorial complex")
      .def("Closeness", &CombinatorialComplex::Closeness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for the given node via nodes "
           "with rank=q")
      .def("ClosenessAll", &CombinatorialComplex::ClosenessAll, "k"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the closeness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("Betweenness", &CombinatorialComplex::Betweenness, "node"_a, "q"_a,
           "weighted"_a = false,
           "Calculate the betweenness centrality for the given node via nodes "
           "with rank=q")
      .def("BetweennessAll", &CombinatorialComplex::BetweennessAll, "k"_a,
           "q"_a, "weighted"_a = false,
           "Calculate the betweenness centrality for all nodes with rank=k via "
           "nodes with rank=q")
      .def("GetEdges", &CombinatorialComplex::GetAll,
           "Get all edges in the combinatorial complex")
      .def("GetElementsWithRank", &CombinatorialComplex::GetElementsWithRank,
           "rank"_a,
           "Get all elements of the combinatorial complex with the specified "
           "rank")
      .def("GetAll", &CombinatorialComplex::GetAll,
           "Get all elements in the combinatorial complex")
      .def("TotalCount", &CombinatorialComplex::TotalCount,
           "Get the total count of elements in the combinatorial complex")
      .def("FVector", &CombinatorialComplex::FVector,
           "Get the F-vector of the combinatorial complex")
      .def("Dimension", &CombinatorialComplex::Dimension,
           "Get the dimension of the combinatorial complex")
      .def("EulerCharacteristic", &CombinatorialComplex::EulerCharacteristic,
           "Get the Euler characteristic of the combinatorial complex")

      .def("BuildFromBinary", &CombinatorialComplex::BuildFromBinary,
           "binary"_a, "on_column"_a = false,
           "Build combinatorial complex from a binary matrix with flag "
           "on_column. "
           "If on_column=true then edges is on column, otherwise on rows")

      .def("Clear", &CombinatorialComplex::Clear,
           "Clear all elements in the combinatorial complex");
}
