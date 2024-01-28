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

void SimplicialComplex::AddArc(const std::vector<int> &from, const std::vector<int> &to) {
    hasse_.AddArc(from, to);
}

void AddCofaces(const std::vector<std::vector<int>> &g, int depth,
                int max_depth, std::vector<int> cur_node, std::vector<int> neighbors,
                SimplicialComplex &simpl) {
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
        simpl.AddArc(cur_node, new_node);
        AddCofaces(g, depth + 1, max_depth, new_node, new_neighbors, simpl);
    }
}

SimplicialComplex SimplicialComplex::CreateCliqueGraph(const std::vector<std::vector<int>> &g, int k, int method) {
    int n = g.size();
    if (n == 0) {
        throw std::runtime_error("Empty graph");
    }
    if (n != g[0].size()) {
        throw std::runtime_error("Wrong format of adjacency matrix");
    }

    SimplicialComplex result;

    if (method == 0) {  // incremental
        for (int v = 0; v < n; v++) {
            std::vector<int> N;
            for (int u = 0; u < v; u++) {
                if (g[v][u]) {
                    N.push_back(u);
                }
            }
            AddCofaces(g, 1, k, {v}, N, result);
        }
    } else {
        for (int v = 0; v < n; v++) {
            for (int u = v + 1; u < n; u++) {
                if (g[v][u]) {
                    result.AddComplex({v, u});
                }
            }
        }

        result.hasse_.DebugPrintAll();
        for (int i = 1; i + 1 <= k; i++) {
        }
        throw std::runtime_error("Not implemented");
    }
    return result;
}

std::vector<std::vector<int>> SimplicialComplex::GetMaxFaces() {
    auto data = hasse_.GetMaxFaces();
    std::vector<std::vector<int>> result;
    for (size_t i = 0; i < data.size(); i++) {
        result.push_back(data[i]->data);
        sort(result.back().begin(), result.back().end());
    }
    return result;
}
