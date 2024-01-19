# SimplicialExtended
Library for simplex complex and etc. Now implemented Hassse diagram.



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