#include <array>
#include <iostream>

#include "SimplicialComplex.h"
#include "generator.h"

// Print only max cliques
void PrintCliques(SimplicialComplex &simpl) {
    auto cliques = simpl.GetMaxFaces();
    for (auto max_cliq : cliques) {
        std::cerr << "{";
        bool is_start = true;
        for (int v : max_cliq) {
            if (!is_start) {
                std::cerr << ", ";
            } else {
                is_start = false;
            }
            std::cerr << v;
        }
        std::cerr << "}\n";
    }
}

void MyTest() {
    int n = 5;
    std::vector<std::vector<int>> edges = {{1, 2},
                                           {2, 3},
                                           {1, 3},
                                           {2, 4},
                                           {4, 5}};
    std::vector<std::vector<int>> g(n, std::vector<int>(n));
    for (const auto &e : edges) {
        int v = e[0] - 1;
        int u = e[1] - 1;
        g[v][u] = g[u][v] = 1;
    }

    auto simpl = std::unique_ptr<SimplicialComplex>(SimplicialComplex::CreateCliqueGraph(g, 5));
    PrintCliques(*simpl);
}

void RandomTest(int n, double prob) {
    auto g = GenerateRandomMatrix(n, prob);

    auto start = std::chrono::high_resolution_clock::now();

    auto simpl = std::unique_ptr<SimplicialComplex>(SimplicialComplex::CreateCliqueGraph(g, 5));

    auto end = std::chrono::high_resolution_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cerr << "Elapsed time of Random test in seconds : " << 1.0 * elapsed / 1000 / 1000 << "\n";
}

int main() {
    // MyTest();
    RandomTest(150, 0.5);
}