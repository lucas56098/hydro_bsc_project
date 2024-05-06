cmake_minimum_required(VERSION 3.5)
project(hydro_bsc_project VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(hydro main.cpp Point.cpp Cell.cpp Mesh.cpp vmp/Halfplane.cpp vmp/VoronoiCell.cpp vmp/VoronoiMesh.cpp)

# Set the name of the compiled program to "hydro"
set_target_properties(hydro PROPERTIES OUTPUT_NAME "hydro")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
