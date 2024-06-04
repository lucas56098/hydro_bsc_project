cmake_minimum_required(VERSION 3.5)
project(hydro_bsc_project VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(hydro 
    
    # main program
    main.cpp 
    Mesh.tpp
    Solver.tpp

    # different cell types
    cell_types/Point.cpp 
    cell_types/Cell.cpp 
    cell_types/Q_Cell.cpp
    cell_types/Conway_Cell.cpp
    cell_types/SWE_Cell.cpp

    # further utilities
    utilities/Functions.cpp
    
    # imported voronoi_mesh_project
    vmp/Halfplane.cpp 
    vmp/VoronoiCell.cpp 
    vmp/VoronoiMesh.cpp
    )

# Set the name of the compiled program to "hydro"
set_target_properties(hydro PROPERTIES OUTPUT_NAME "hydro")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
