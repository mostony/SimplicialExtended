cmake_minimum_required(VERSION 3.11...3.26)
project(simpl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

add_subdirectory(testing)

Include(FetchContent)

FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.11.1
)
FetchContent_MakeAvailable(pybind11)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
)

FetchContent_MakeAvailable(Catch2)

add_executable(tests test.cpp src/Hasse.cpp
        src/HyperGraph.cpp
        src/Node.cpp
        src/SimplicialComplex.cpp
        src/CombinatorialComplex.cpp)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

pybind11_add_module(simpl
        module_creation.cpp
        src/Hasse.cpp
        src/HyperGraph.cpp
        src/Node.cpp
        src/SimplicialComplex.cpp
        src/CombinatorialComplex.cpp)
