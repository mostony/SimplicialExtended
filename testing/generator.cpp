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

#include <iostream>

#include <filesystem>
#include <fstream>

std::vector<std::vector<int>> GenerateImageDatasetMatrix(
    int n, int k, int random_seed = 228) {
    // TODO
    // auto generator = py::module::import("generate_image_dataset");
    // auto resultobj = generator.attr("GenerateImageDataset")(n, k, random_seed);
    // return resultobj.cast<std::vector<std::vector<int>>>();
    std::string call_string = "python3 ./testing/generate_image_dataset.py " +
                              std::to_string(n) + " " + std::to_string(k) +
                              " " + std::to_string(random_seed);

    std::cerr << "calling: " << call_string << "\n";
    system(call_string.c_str());
    std::cerr << std::filesystem::current_path() << "\n";

    std::ifstream in("testing/images.txt");  // TODO very bad
    std::vector<std::vector<int>> g(n, std::vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int x;
            in >> x;
            g[i][j] = x;
        }
    }
    return g;
}
