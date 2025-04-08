cmake_minimum_required(VERSION 3.11...3.26)
project(simpl)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

add_subdirectory(testing)

Include(FetchContent)

FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.0
)

FetchContent_MakeAvailable(Catch2)

FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
)

FetchContent_MakeAvailable(Eigen3)

FetchContent_Declare(
        Spectra
        GIT_REPOSITORY https://github.com/yixuan/spectra.git
        GIT_TAG v1.1.0
)

FetchContent_MakeAvailable(Spectra)

add_executable(tests test.cpp src/Hasse.cpp
        src/HyperGraph.cpp
        src/Node.cpp
        src/SimplicialComplex.cpp
        src/CombinatorialComplex.cpp
        src/Graph.cpp
        src/AbstractModel.cpp
)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(tests PRIVATE Eigen3::Eigen)
target_link_libraries(tests PRIVATE Spectra)

pybind11_add_module(simpl
        module_creation.cpp
        src/Hasse.cpp
        src/HyperGraph.cpp
        src/Node.cpp
        src/SimplicialComplex.cpp
        src/CombinatorialComplex.cpp
        src/Graph.cpp
        src/AbstractModel.cpp
)

target_link_libraries(simpl PRIVATE Eigen3::Eigen)
target_link_libraries(simpl PRIVATE Spectra)

find_package(OpenMP REQUIRED)

target_link_libraries(simpl PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(tests PUBLIC OpenMP::OpenMP_CXX)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(simpl PUBLIC ${ARMADILLO_LIBRARIES})

target_link_libraries(tests PUBLIC ${ARMADILLO_LIBRARIES})