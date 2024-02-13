cmake_minimum_required(VERSION 3.5...3.26)
project(simpl)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

# add_compile_options(-fsanitize=thread)
# add_link_options(-fsanitize=thread)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(extern/pybind11)

# pybind11_add_module(simpl 
#                 module_creation.cpp 
#                 Hasse.cpp
#                 HyperGraph.cpp
#                 Node.cpp
#                 SimplicialComplex.cpp)

add_subdirectory(src)
add_subdirectory(testing)


add_executable(simplex_test test_clique_matrix.cpp
        src/Hasse.cpp src/Hasse.h
        src/SimplicialComplex.cpp src/SimplicialComplex.h                    
        src/Node.cpp src/Node.h
        testing/generator.cpp testing/generator.h
)

# target_link_libraries(simplex_test
#     pybind11::embed
# )


