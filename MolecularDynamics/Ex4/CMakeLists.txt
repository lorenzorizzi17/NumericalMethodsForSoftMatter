cmake_minimum_required(VERSION 3.10)

project(NDimMolDyn CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -fopenmp")

find_package(SFML 2.5 COMPONENTS system window graphics QUIET)
find_package(OpenMP REQUIRED)

add_executable(MultipleChains.out
    code/multipleMarkovChains.cpp
)
target_link_libraries(MultipleChains.out
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_executable(2DimMolDyn.out
    code/main.cpp
)
target_link_libraries(2DimMolDyn.out
    PRIVATE
)
