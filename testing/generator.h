#pragma once

#include <chrono>
#include <random>
#include <vector>

// #include "pybind11/embed.h"

// namespace py = pybind11;

std::vector<std::vector<int>> GenerateRandomMatrix(int n, double p);

std::vector<std::vector<int>> GenerateImageDatasetMatrix(int n, int k,
                                                         int random_seed);
