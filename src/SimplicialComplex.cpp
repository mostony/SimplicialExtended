#include "SimplicialComplex.h"

void SimplicialComplex::AddComplex(std::vector<VertexId> complex) {
    hasse_.RecursiveAddNode(complex);
}

void SimplicialComplex::RemoveComplex(std::vector<VertexId> complex) {
    hasse_.RemoveNode(complex);
}

void SimplicialComplex::Debug() {
    hasse_.DebugPrintAll();
}

int SimplicialComplex::HasseSize() {
    return hasse_.Size();
}

void SimplicialComplex::AddArc(const std::vector<int>& from,
                               const std::vector<int>& to) {
    hasse_.AddArc(from, to);
}

void AddCofaces(const std::vector<std::vector<int>>& g, int depth,
                int max_depth, std::vector<int> cur_node,
                std::vector<int> neighbors, SimplicialComplex* simpl) {
    if (depth > max_depth) {
        return;
    }
    // simpl.AddComplex(cur_node);
    int n = g.size();
    for (int u : neighbors) {
        int ptr = 0;
        std::vector<int> new_neighbors;
        std::vector<int> new_node = cur_node;
        new_node.push_back(u);

        // intersection of 2 vectors
        for (int x = 0; x < u; x++) {
            if (g[u][x]) {
                while (ptr < neighbors.size() && neighbors[ptr] < x) {
                    ptr += 1;
                }
                if (ptr == neighbors.size()) {
                    break;
                }
                if (neighbors[ptr] == x) {
                    new_neighbors.push_back(x);
                }
            }
        }
        simpl->AddArc(cur_node, new_node);
        AddCofaces(g, depth + 1, max_depth, new_node, new_neighbors, simpl);
    }
}

SimplicialComplex* SimplicialComplex::CreateCliqueGraph(
    const std::vector<std::vector<int>>& g, int k, int method,
    int total_threads) {
    int n = g.size();
    if (n == 0) {
        throw std::runtime_error("Empty graph");
    }
    if (n != g[0].size()) {
        throw std::runtime_error("Wrong format of adjacency matrix");
    }

    auto result = new SimplicialComplex();
    if (total_threads == -1) {
        // TODO I'm not sure about this line
        total_threads = std::thread::hardware_concurrency();
    }

    if (total_threads == 0) {
        throw std::runtime_error("Zero threads");
    }

    if (method == 0) {  // incremental
        std::vector<std::thread> threads;
        std::vector<SimplicialComplex*> thread_results;

        auto AddSubsetVertices = [&](int l, int r,
                                     SimplicialComplex* thread_result) {
            for (int i = l; i < r; i++) {
                int v = i;
                std::vector<int> N;
                for (int u = 0; u < v; u++) {
                    if (g[v][u]) {
                        N.push_back(u);
                    }
                }
                std::vector<int> cur_node = {v};
                AddCofaces(g, 1, k, cur_node, N, thread_result);
            }
        };

        int bucket = std::max((n + total_threads - 1) / total_threads, 1);
        for (int index_thread = 0; index_thread < total_threads;
             index_thread++) {
            int l = index_thread * bucket;
            int r = std::min((index_thread + 1) * bucket, n);
            // std::cerr << l << ' ' << r << std::endl;
            if (l < r) {
                thread_results.push_back(new SimplicialComplex());
                std::thread th(AddSubsetVertices, l, r, thread_results.back());
                threads.emplace_back(std::move(th));
            }
        }

        for (auto& thread : threads) {
            thread.join();
        }

        for (auto thread_result : thread_results) {
            Merge(result, thread_result);
        }
    } else {
        std::vector<std::thread> threads;
        std::vector<SimplicialComplex*> thread_results;

        auto AddSubsetVertices = [&](int l, int r,
                                     SimplicialComplex* thread_result) {
            std::vector<std::vector<int>> cur_layer;
            for (int v = l; v < r; v++) {
                std::vector<int> from = {v};
                for (int u = 0; u < v; u++) {
                    if (g[v][u]) {
                        std::vector<int> to = {v, u};
                        thread_result->AddArc(from, to);
                    }
                }
                cur_layer.push_back(from);
            }
            for (int it = 1; it + 1 < k; it++) {
                std::vector<std::vector<int>> next_layer;

                for (auto middle : cur_layer) {
                    auto middle_node = thread_result->hasse_.GetNode(middle);
                    for (auto sigma1 : middle_node->sons) {
                        int v1 = sigma1.back();
                        for (auto sigma2 : middle_node->sons) {
                            int v2 = sigma2.back();
                            if (v1 > v2 && g[v1][v2]) {
                                auto new_data = sigma1;
                                new_data.push_back(v2);
                                thread_result->AddArc(sigma1, new_data);
                            }
                        }
                        next_layer.push_back(sigma1);
                    }
                }
                std::swap(cur_layer, next_layer);
            }
        };

        int bucket = std::max((n + total_threads - 1) / total_threads, 1);
        for (int index_thread = 0; index_thread < total_threads;
             index_thread++) {
            int l = index_thread * bucket;
            int r = std::min((index_thread + 1) * bucket, n);
            if (l < r) {
                thread_results.push_back(new SimplicialComplex());
                std::thread th(AddSubsetVertices, l, r, thread_results.back());
                threads.emplace_back(std::move(th));
            }
        }

        for (auto& thread : threads) {
            thread.join();
        }

        for (auto thread_result : thread_results) {
            Merge(result, thread_result);
        }
    }
    return result;
}

std::vector<std::vector<int>> SimplicialComplex::GetMaxFaces() {
    auto data = hasse_.GetMaxFaces();
    std::vector<std::vector<int>> result;
    for (size_t i = 0; i < data.size(); i++) {
        result.push_back(data[i]->data);
        std::sort(result.back().begin(), result.back().end());
    }
    return result;
}

void Merge(SimplicialComplex* current, SimplicialComplex* other) {
    Merge(current->hasse_, other->hasse_);
    delete other;
}
