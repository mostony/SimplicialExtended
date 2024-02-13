#include <fstream>
#include "generator.h"

const int TRY = 1;

int main() {
    for (int k : {5, 10, 15, 20}) {
        for (int n = 32; n <= 32; n *= 2) {
            for (int random_seed = 0; random_seed < TRY; random_seed++) {
                GenerateImageDatasetMatrix(std::min(n, 20000), k, random_seed);
            }
            // std::ofstream fout("./testing/times" + std::to_string(k) + ".txt",
            //                    std::ios_base::app);
            // fout << "\n";
        }
    }
}