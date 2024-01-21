#include "generator.h"

// fixed state
std::mt19937 rng(5151);
std::uniform_real_distribution<double> uniform(0, 1);

// generate graph on n vertex
// each edge has probability of exists equals p
std::vector<std::vector<int>> GenerateRandomMatrix(int n, double p) {
    std::vector<std::vector<int>> g(n, std::vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (uniform(rng) < p) {
                g[i][j] = g[j][i] = 1;
            }
        }
    }
    return g;
}