cmake_minimum_required(VERSION 3.26)
project(simpl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add_subdirectory(external/pybind11)
# pybind11_add_module(simpl 
#                 module_creation.cpp 
#                 Hasse.cpp
#                 HyperGraph.cpp
#                 Node.cpp
#                 SimplicialComplex.cpp)

add_executable(simplex_test test_clique_matrix.cpp
        Hasse.cpp Hasse.h
        SimplicialComplex.cpp SimplicialComplex.h                    
        Node.cpp Node.h
        generator.cpp generator.h
)


