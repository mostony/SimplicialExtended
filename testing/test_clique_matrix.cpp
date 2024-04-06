#include <fstream>
#include <iostream>

#include "../src/SimplicialComplex.h"
#include "generator.h"

// Print only max cliques
void PrintCliques(SimplicialComplex& simpl) {
  auto cliques = simpl.GetMaxSimplices();
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
  std::vector<std::vector<int>> edges = {
      {1, 2}, {2, 3}, {1, 3}, {2, 4}, {4, 5}};
  std::vector<std::vector<int>> g(n, std::vector<int>(n));
  for (const auto& e : edges) {
    int v = e[0] - 1;
    int u = e[1] - 1;
    g[v][u] = g[u][v] = 1;
  }

  auto simpl = std::unique_ptr<SimplicialComplex>(
      SimplicialComplex::CreateCliqueGraph(g, 5, 1));
  PrintCliques(*simpl);
}

void RandomTest(int n, double prob) {
  auto g = GenerateRandomMatrix(n, prob);

  auto start = std::chrono::high_resolution_clock::now();

  auto simpl = std::unique_ptr<SimplicialComplex>(
      SimplicialComplex::CreateCliqueGraph(g, 5, 0, -1));

  auto end = std::chrono::high_resolution_clock::now();

  auto elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();
  std::cerr << "Elapsed time of Random test in seconds : "
            << 1.0 * elapsed / 1000 / 1000 << "\n";
}

void LetterTest5() {
  auto g = GenerateImageDatasetMatrix(20000, 5, 1);
  for (int method : {0, 1}) {
    for (int threads : {-1, 1}) {
      auto start = std::chrono::high_resolution_clock::now();

      auto simpl = std::unique_ptr<SimplicialComplex>(
          SimplicialComplex::CreateCliqueGraph(g, 50, method, threads));

      auto end = std::chrono::high_resolution_clock::now();

      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      if (method == 0) {
        std::cerr << "Incremental method\n";
      } else {
        std::cerr << "Rieser method\n";
      }
      if (threads == -1) {
        std::cerr << "All threads\n";
      } else {
        std::cerr << "Single thread\n";
      }
      std::cerr << "Elapsed time of letter_test for "
                << "k = 5"
                << " in seconds : " << 1.0 * elapsed / 1000 / 1000 << "\n";
      std::cerr << "Hasse size is : " << simpl->HasseSize() << "\n";
    }
  }
}

void LetterTest10() {
  auto g = GenerateImageDatasetMatrix(20000, 10, 0);
  for (int method : {0, 1}) {
    for (int threads : {-1, 1}) {
      auto start = std::chrono::high_resolution_clock::now();

      auto simpl = std::unique_ptr<SimplicialComplex>(
          SimplicialComplex::CreateCliqueGraph(g, 50, method, threads));

      auto end = std::chrono::high_resolution_clock::now();

      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      if (method == 0) {
        std::cerr << "Incremental method\n";
      } else {
        std::cerr << "Rieser method\n";
      }
      if (threads == -1) {
        std::cerr << "All threads\n";
      } else {
        std::cerr << "Single thread\n";
      }
      std::cerr << "Elapsed time of letter_test for "
                << "k = 10"
                << " in seconds : " << 1.0 * elapsed / 1000 / 1000 << "\n";
      std::cerr << "Hasse size is : " << simpl->HasseSize() << "\n";
    }
  }
}

void LetterTest15() {
  auto g = GenerateImageDatasetMatrix(20000, 15, 0);
  for (int method : {0, 1}) {
    for (int threads : {-1, 1}) {
      auto start = std::chrono::high_resolution_clock::now();

      auto simpl = std::unique_ptr<SimplicialComplex>(
          SimplicialComplex::CreateCliqueGraph(g, 50, method, threads));

      auto end = std::chrono::high_resolution_clock::now();

      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      if (method == 0) {
        std::cerr << "Incremental method\n";
      } else {
        std::cerr << "Rieser method\n";
      }
      if (threads == -1) {
        std::cerr << "All threads\n";
      } else {
        std::cerr << "Single thread\n";
      }
      std::cerr << "Elapsed time of letter_test for "
                << "k = 15"
                << " in seconds : " << 1.0 * elapsed / 1000 / 1000 << "\n";
      std::cerr << "Hasse size is : " << simpl->HasseSize() << "\n";
    }
  }
}

void LetterTest20() {
  auto g = GenerateImageDatasetMatrix(20'000, 20, 0);

  for (int method : {0, 1}) {
    for (int threads : {-1, +1}) {
      auto start = std::chrono::high_resolution_clock::now();

      auto simpl = std::unique_ptr<SimplicialComplex>(
          SimplicialComplex::CreateCliqueGraph(g, 50, method, threads));

      auto end = std::chrono::high_resolution_clock::now();

      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      if (method == 0) {
        std::cerr << "Incremental method\n";
      } else {
        std::cerr << "Rieser method\n";
      }
      if (threads == -1) {
        std::cerr << "All threads\n";
      } else {
        std::cerr << "Single thread\n";
      }
      std::cerr << "Elapsed time of letter_test for "
                << "k = 20"
                << " in seconds : " << 1.0 * elapsed / 1000 / 1000 << "\n";
      std::cerr << "Hasse size is : " << simpl->HasseSize() << "\n";
    }
  }
}

const int TRY = 1;

void Test() {
  for (int method : {0, 1}) {
    for (int threads : {-1, +1}) {
      for (int k : {5, 10, 15, 20}) {
        std::string path = "./testing/algo" + std::to_string(method) +
                           std::to_string(threads == -1) + std::to_string(k) +
                           ".txt";
        for (int n = 32; n <= 2048; n *= 2) {
          double mean_time = 0;
          for (int random_seed = 0; random_seed < TRY; random_seed++) {
            auto g =
                GenerateImageDatasetMatrix(std::min(n, 20'000), k, random_seed);

            auto start = std::chrono::high_resolution_clock::now();
            auto simpl = std::unique_ptr<SimplicialComplex>(
                SimplicialComplex::CreateCliqueGraph(g, 50, method, threads));

            auto end = std::chrono::high_resolution_clock::now();

            auto elapsed =
                std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                      start)
                    .count();

            mean_time += 1.0 * elapsed / 1000 / 1000;
          }
          mean_time /= TRY;
          std::cerr << "Mean Elapsed time " << method << " " << threads << " "
                    << mean_time << "\n";

          std::ofstream fout(path, std::ios_base::app);
          fout << mean_time << " ";
        }
      }
    }
  }
}

int main() {
  Test();
  // MyTest();
  // RandomTest(120, 0.5);
  // LetterTest5();
  // LetterTest10();
  // LetterTest15();
  // LetterTest20();
}