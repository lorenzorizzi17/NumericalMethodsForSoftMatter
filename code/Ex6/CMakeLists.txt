cmake_minimum_required(VERSION 3.10)

project(NDimMolDyn CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

find_package(SFML 2.5 COMPONENTS system window graphics REQUIRED)

#add_executable(2DimMolDyn.out
#    code/sfml.cpp
#)
#target_link_libraries(2DimMolDyn.out
#    PRIVATE
#    sfml-system
#    sfml-window
#    sfml-graphics
#)

add_executable(3DimMolDyn.out
    code/main.cpp
)
target_link_libraries(3DimMolDyn.out
    PRIVATE
)
