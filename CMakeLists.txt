cmake_minimum_required(VERSION 3.10)
project(paradiso)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

include_directories(
        /opt/intel/oneapi/mkl/latest/include
        include
)

link_directories(
        /opt/intel/oneapi/mkl/latest/lib/intel64
)

add_executable(solver
        paradiso_solver.cpp
        include/MatrixLoader.cpp
)

target_link_libraries(solver
        mkl_intel_ilp64
        mkl_intel_thread
        mkl_core
        iomp5
        pthread
        m
        dl
        OpenMP::OpenMP_CXX
)

add_executable(generator
        matrixGenerator.cpp
)