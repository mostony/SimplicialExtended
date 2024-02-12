cmake_minimum_required(VERSION 3.5...3.26)

# add_executable(test_generating test_generating.cpp generator.h generator.cpp)


# Copy python file to builde project
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/generate_image_dataset.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/testing)
