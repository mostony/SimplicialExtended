#include <fstream>
#include <iostream>
#include <string>

int main() {
  for (int k : {5, 10, 15, 20}) {
    std::ifstream in("times" + std::to_string(k) + ".txt");
    std::ofstream out("mean_times" + std::to_string(k) + ".txt");

    for (int n = 32; n <= 32 * 1024; n *= 2) {
      double mean_time = 0;
      for (int it = 0; it < 10; it++) {
        double x;
        in >> x;
        mean_time += x;
      }
      mean_time /= 10;
      out << mean_time << ", ";
    }
  }
}