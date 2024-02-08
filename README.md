# SimplicialExtended
Library for simplex complex and etc. Now implemented using Hasse diagram.



## Usage

### Installing pybind11
```
git submodule add -b stable ../../pybind/pybind11 external/pybind11
git submodule update --init
```

### Create module
```
cmake -S . -B build
cd build
cmake --build .
```

then just import simpl.

### Testing
For testing you can generate letter_image dataset using /dataset/generate_image_dataset.py.
Then launch test_clique_matrix.cpp