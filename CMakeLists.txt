cmake_minimum_required(VERSION 3.5...3.26)
project(simpl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(extern/pybind11)

add_subdirectory(src)

add_subdirectory(testing)

pybind11_add_module(simpl
        module_creation.cpp
        src/Hasse.cpp
        src/HyperGraph.cpp
        src/Node.cpp
        src/SimplicialComplex.cpp)
