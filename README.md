# SimplicialExtended
Library for creating complex networks and calculating features on them.

Supported networks: simplicial complexes, graphs, hypergraphs and combinatorial complexes

Main supported features: Betti number, custom function addition and removing, centralities, feature matrix, laplacian eigen values.

## Usage

```
cmake -S . -B build
cd build
cmake --build .
```

This generate unit test and python module (simpl). This module can be latter imported to python code.

Example of usage [here](tutorial/example.ipynb)