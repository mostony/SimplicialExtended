#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "src/CombinatorialComplex.h"
#include "src/Graph.h"
#include "src/SimplicialComplex.h"

TEST_CASE("Simplicial Complex incidence") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3, 4});
  simpl.AddComplex({1, 5, 6});

  // v = 1

  // k = 2
  {
    auto incidence = simpl.Incidence({1}, 2);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1, 4}, {1, 3}, {1, 2}, {1, 6}, {1, 5}})));
  }

  // k = 3
  {
    auto incidence = simpl.Incidence({1}, 3);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1, 3, 4}, {1, 2, 4}, {1, 2, 3}, {1, 5, 6}})));
  }

  // k = 4
  {
    auto incidence = simpl.Incidence({1}, 4);
    REQUIRE_THAT(incidence, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{1, 2, 3, 4}})));
  }
}

TEST_CASE("Simplicial Complex degree") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3, 4});
  simpl.AddComplex({1, 5, 6});

  // k = 2
  {
    auto degree = simpl.Degree({1}, 2);
    REQUIRE_THAT(degree,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1}, {2}, {3}, {4}, {5}, {6}})));
  }
  // k = 3
  {
    auto degree = simpl.Degree({1}, 3);
    REQUIRE_THAT(degree,
                 Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                     {{1}, {2}, {3}, {4}, {5}, {6}})));
  }

  // k = 4
  {
    auto degree = simpl.Degree({1}, 4);
    sort(degree.begin(), degree.end());
    REQUIRE_THAT(degree,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{1}, {2}, {3}, {4}})));
  }
}

TEST_CASE("Combinatorial Complex incidence") {
  CombinatorialComplex comb;
  comb.Build(
      {{1}, {2}, {3}, {4}, {1, 2}, {1, 4}, {2, 4}, {1, 2, 3}, {2, 3, 4}});
  // v = 1, k = 2
  {
    auto incidence = comb.Incidence({1}, 2);
    REQUIRE_THAT(incidence,
                 Catch::Matchers::UnorderedEquals(
                     std::vector<std::vector<int>>({{1, 2}, {1, 4}})));
  }

  // v = 1, k = 3
  {
    auto incidence = comb.Incidence({1}, 3);
    REQUIRE_THAT(incidence, Catch::Matchers::UnorderedEquals(
                                std::vector<std::vector<int>>({{1, 2, 3}})));
  }
}

TEST_CASE("Combinatorial Complex degree") {
  CombinatorialComplex comb;
  comb.Build(
      {{1}, {2}, {3}, {4}, {1, 2}, {1, 4}, {2, 4}, {1, 2, 3}, {2, 3, 4}});

  // v = 1, k = 2
  {
    auto degree = comb.Degree({1}, 2);
    REQUIRE_THAT(degree, Catch::Matchers::UnorderedEquals(
                             std::vector<std::vector<int>>({{1}, {2}, {4}})));
  }
  // v = 1, k = 3
  {
    auto degree = comb.Degree({1}, 3);
    REQUIRE_THAT(degree, Catch::Matchers::UnorderedEquals(
                             std::vector<std::vector<int>>({{1}, {2}, {3}})));
  }
}

TEST_CASE("Betti number") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3});
  simpl.AddComplex({3, 4});
  simpl.AddComplex({4, 5, 6});
  simpl.AddComplex({4, 5, 7});
  simpl.AddComplex({6, 7});
  simpl.AddComplex({2, 7});
  simpl.RemoveComplex({1, 2});

  REQUIRE(simpl.BettiNumber(0) == 1);
  REQUIRE(simpl.BettiNumber(1) == 2);
  REQUIRE(simpl.BettiNumber(2) == 0);
}

TEST_CASE("Simplicial closeness") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3, 4});
  simpl.AddComplex({1, 5, 6});

  // max_rank = 2
  REQUIRE(simpl.Closeness({1}, 2) == Catch::Approx(1.0));
  REQUIRE(simpl.Closeness({2}, 2) == Catch::Approx(0.71428571428));
  REQUIRE(simpl.Closeness({5}, 2) == Catch::Approx(0.625));

  REQUIRE(simpl.Closeness({1}, 3) == Catch::Approx(1));
  REQUIRE(simpl.Closeness({1}, 4) == Catch::Approx(0));
}

TEST_CASE("Simplicial Betweenness") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3, 4});
  simpl.AddComplex({1, 5, 6});

  // max_rank = 2
  REQUIRE(simpl.Betweenness({1}, 2) == Catch::Approx(0.6));
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

  // max_rank = 2
  REQUIRE(graph.Betweenness({1}, 2) == Catch::Approx(0.65));
}

TEST_CASE("KEKW") {
  SimplicialComplex simpl;
  simpl.AddComplex({1, 2, 3});
  simpl.AddComplex({3, 4});
  simpl.AddComplex({4, 5, 6});
  simpl.AddComplex({4, 5, 7});
  simpl.AddComplex({6, 7});
  simpl.AddComplex({2, 7});
  simpl.RemoveComplex({1, 2});
  simpl.GetMaxSimplices();
}
