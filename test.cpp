#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "src/CombinatorialComplex.h"
#include "src/Graph.h"
#include "src/Hasse.h"
#include "src/HyperGraph.h"
#include "src/SimplicialComplex.h"
#include <Eigen/Dense>

TEST_CASE("Simplicial Complex incidence") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  // v = 1, k = 1
  {
    auto incidence = simpl.Incidence({1}, 1);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1, 4}, {1, 3}, {1, 2}, {1, 6}, {1, 5}})));
    auto incidence_degree = simpl.IncidenceDegree({1}, 1);
    REQUIRE(incidence_degree == 5);
  }

  // v = 1, k = 2
  {
    auto incidence = simpl.Incidence({1}, 2);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1, 3, 4}, {1, 2, 4}, {1, 2, 3}, {1, 5, 6}})));
    auto incidence_degree = simpl.IncidenceDegree({1}, 2);
    REQUIRE(incidence_degree == 4);
  }

  // v = 1, k = 3
  {
    auto incidence = simpl.Incidence({1}, 3);
    REQUIRE_THAT(incidence, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{1, 2, 3, 4}})));
    auto incidence_degree = simpl.IncidenceDegree({1}, 3);
    REQUIRE(incidence_degree == 1);
  }
}

TEST_CASE("Simplicial Complex degree") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  // v = 1, k = 1
  {
    auto adjacency = simpl.Adjacency({1}, 1);
    REQUIRE_THAT(adjacency,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{2}, {3}, {4}, {5}, {6}})));
    auto degree = simpl.Degree({1}, 1);
    REQUIRE(degree == 5);
  }

  // v = 1, k = 2
  {
    auto adjacency = simpl.Adjacency({1}, 2);
    REQUIRE_THAT(adjacency,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{2}, {3}, {4}, {5}, {6}})));
    auto degree = simpl.Degree({1}, 2);
    REQUIRE(degree == 5);
  }

  // v = 1, k = 3
  {
    auto adjacency = simpl.Adjacency({1}, 3);
    REQUIRE_THAT(adjacency,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{2}, {3}, {4}})));
    auto degree = simpl.Degree({1}, 3);
    REQUIRE(degree == 3);
  }
}

TEST_CASE("Combinatorial Complex incidence") {
  CombinatorialComplex comb;
  comb.Build(
      {{1}, {2}, {3}, {4}, {1, 2}, {1, 4}, {2, 4}, {1, 2, 3}, {2, 3, 4}});

  // v = 1, k = 1
  {
    auto incidence = comb.Incidence({1}, 1);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{1, 2}, {1, 4}})));
    auto incidence_degree = comb.IncidenceDegree({1}, 1);
    REQUIRE(incidence_degree == 2);
  }

  // v = 1, k = 2
  {
    auto incidence = comb.Incidence({1}, 2);
    REQUIRE_THAT(incidence, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{1, 2, 3}})));
    auto incidence_degree = comb.IncidenceDegree({1}, 2);
    REQUIRE(incidence_degree == 1);
  }
}

TEST_CASE("Combinatorial Complex degree") {
  CombinatorialComplex comb;
  comb.Build(
      {{1}, {2}, {3}, {4}, {1, 2}, {1, 4}, {2, 4}, {1, 2, 3}, {2, 3, 4}});

  // v = 1, k = 1
  {
    auto adjacency = comb.Adjacency({1}, 1);
    REQUIRE_THAT(adjacency, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{2}, {4}})));
    auto degree = comb.Degree({1}, 1);
    REQUIRE(degree == 2);
  }

  // v = 1, k = 2
  {
    auto adjacency = comb.Adjacency({1}, 2);
    REQUIRE_THAT(adjacency, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{2}, {3}})));
    auto degree = comb.Degree({1}, 2);
    REQUIRE(degree == 2);
  }
}

TEST_CASE("Betti numbers") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3});
  simpl.AddSimplex({3, 4});
  simpl.AddSimplex({4, 5, 6});
  simpl.AddSimplex({4, 5, 7});
  simpl.AddSimplex({6, 7});
  simpl.AddSimplex({2, 7});
  simpl.RemoveSimplex({1, 2});

  REQUIRE(simpl.BettiNumber(0) == 1);
  REQUIRE(simpl.BettiNumber(1) == 2);
  REQUIRE(simpl.BettiNumber(2) == 0);
}

TEST_CASE("Simplicial closeness") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  REQUIRE(simpl.Closeness({1}, 1) == Catch::Approx(1.0));
  REQUIRE(simpl.Closeness({2}, 1) == Catch::Approx(0.71428571428));
  REQUIRE(simpl.Closeness({5}, 1) == Catch::Approx(0.625));

  REQUIRE(simpl.Closeness({1}, 2) == Catch::Approx(1));
  REQUIRE(simpl.Closeness({1}, 3) == Catch::Approx(0.6));

  // weighted = 1.0 check
  REQUIRE(simpl.Closeness({1}, 1, true) == Catch::Approx(1.0));
  REQUIRE(simpl.Closeness({2}, 1, true) == Catch::Approx(0.71428571428));
  REQUIRE(simpl.Closeness({5}, 1, true) == Catch::Approx(0.625));

  REQUIRE(simpl.Closeness({1}, 2, true) == Catch::Approx(1));
  REQUIRE(simpl.Closeness({1}, 3, true) == Catch::Approx(0.6));
}

TEST_CASE("Simplicial Betweenness") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  // v = 1, max_rank = 1
  REQUIRE(simpl.Betweenness({1}, 1) == Catch::Approx(0.6));

  // weighted
  REQUIRE(simpl.Betweenness({1}, 1, true) == Catch::Approx(0.6));
}

TEST_CASE("Graph Betweenness") {
  Graph graph;
  graph.AddEdge(1, 2);
  graph.AddEdge(2, 3);
  graph.AddEdge(3, 4);
  graph.AddEdge(1, 4);
  graph.AddEdge(1, 5);
  graph.AddEdge(1, 6);
  graph.AddEdge(5, 6);

  // v = 1, max_rank = 1
  REQUIRE(graph.Betweenness({1}, 1) == Catch::Approx(0.65));
}

TEST_CASE("Simplex dimensions") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  REQUIRE(simpl.Dimension() == 3);
  simpl.RemoveSimplex({1, 2, 3, 4});
  REQUIRE(simpl.Dimension() == 2);
  simpl.RemoveSimplex({1, 5, 6});
  REQUIRE(simpl.Dimension() == 2);
  simpl.RemoveSimplex({1, 2, 4});
  simpl.RemoveSimplex({1, 2, 3});
  simpl.RemoveSimplex({1, 3, 4});
  simpl.RemoveSimplex({2, 3, 4});
  REQUIRE(simpl.Dimension() == 1);
}

TEST_CASE("Simplex fvector & total count") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});

  REQUIRE_THAT(
      simpl.FVector(),
      Catch::Matchers::UnorderedEquals(
          std::vector<std::pair<int, int>>({{0, 6}, {1, 9}, {2, 5}, {3, 1}})));
  REQUIRE(simpl.TotalCount() == 21);
}

TEST_CASE("Simplex Euler Characteristic") {
  SimplicialComplex simpl;
  simpl.AddSimplex({1, 2, 3});
  simpl.AddSimplex({3, 4});
  simpl.AddSimplex({4, 5, 6});
  simpl.AddSimplex({4, 5, 7});
  simpl.AddSimplex({6, 7});
  simpl.AddSimplex({2, 7});
  simpl.RemoveSimplex({1, 2});

  REQUIRE(simpl.EulerCharacteristic() == -1);
  simpl.Clear();

  simpl.AddSimplex({1, 2, 3, 4});
  simpl.AddSimplex({1, 5, 6});
  REQUIRE(simpl.EulerCharacteristic() == +1);
}

TEST_CASE("Check sign of permutation") {
  REQUIRE(Hasse::CalculateSign({1}, {1, 2, 3, 4}) == 1);
  REQUIRE(Hasse::CalculateSign({2}, {1, 2, 3, 4}) == -1);
  REQUIRE(Hasse::CalculateSign({3}, {1, 2, 3, 4}) == 1);
  REQUIRE(Hasse::CalculateSign({4}, {1, 2, 3, 4}) == -1);
  REQUIRE(Hasse::CalculateSign({2, 4}, {1, 2, 3, 4}) == -1);
  REQUIRE(Hasse::CalculateSign({1, 2, 3}, {1, 2, 3, 4}) == 1);
  REQUIRE(Hasse::CalculateSign({2, 3, 4}, {1, 2, 3, 4}) == -1);

  REQUIRE(Hasse::CalculateSign({2}, {2, 3}) == 1);
  REQUIRE(Hasse::CalculateSign({3}, {2, 3}) == -1);
}
