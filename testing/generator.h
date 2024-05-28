#pragma once

#include <vector>

std::vector<std::vector<int>> GenerateRandomMatrix(int n, double p);

std::vector<std::vector<int>> GenerateImageDatasetMatrix(int n, int k,
                                                         int random_seed);
