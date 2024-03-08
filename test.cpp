#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "src/SimplicialComplex.h"

TEST_CASE("Simplicial incidence test") {
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3, 4});
    simpl.AddComplex({1, 5, 6});

    // k = 2
    {
        auto incidence = simpl.Incidence({1}, 2);
        REQUIRE_THAT(
            incidence,
            Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                {{1, 4}, {1, 3}, {1, 2}, {1, 6}, {1, 5}})));
    }

    // k == 3
    {
        auto incidence = simpl.Incidence({1}, 3);
        REQUIRE_THAT(
            incidence,
            Catch::Matchers::UnorderedEquals(std::vector<std::vector<int>>(
                {{1, 3, 4}, {1, 2, 4}, {1, 2, 3}, {1, 5, 6}})));
    }

    // k == 4
    {
        auto incidence = simpl.Incidence({1}, 4);
        REQUIRE_THAT(incidence,
                     Catch::Matchers::UnorderedEquals(
                         std::vector<std::vector<int>>({{1, 2, 3, 4}})));
    }
}

TEST_CASE("Simplicial degree test") {
    SimplicialComplex simpl;
    simpl.AddComplex({1, 2, 3, 4});
    simpl.AddComplex({1, 5, 6});

    // k == 2
    {
        auto degree = simpl.Degree({1}, 2);
        REQUIRE_THAT(
            degree,
            Catch::Matchers::UnorderedEquals(
                std::vector<std::vector<int>>({{1}, {2}, {3}, {4}, {5}, {6}})));
    }
    // k == 3
    {
        auto degree = simpl.Degree({1}, 3);
        REQUIRE_THAT(
            degree,
            Catch::Matchers::UnorderedEquals(
                std::vector<std::vector<int>>({{1}, {2}, {3}, {4}, {5}, {6}})));
    }

    // k == 4
    {
        auto degree = simpl.Degree({1}, 4);
        sort(degree.begin(), degree.end());
        REQUIRE_THAT(degree,
                     Catch::Matchers::UnorderedEquals(
                         std::vector<std::vector<int>>({{1}, {2}, {3}, {4}})));
    }
}
