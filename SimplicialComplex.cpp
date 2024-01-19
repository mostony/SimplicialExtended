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
