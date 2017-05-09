cmake_minimum_required(VERSION 2.8)

project(internal_branches)

add_compile_options(-O3)
include_directories(ib_eigen)
include_directories(ib_tclap)
include_directories(ib_boost)

add_executable(internal_branches internal_branches.cpp)
target_include_directories(internal_branches PUBLIC ${CMAKE_SOURCE_DIR})

add_executable(FTEC_branchlength FTEC_branchlength.cpp)
target_include_directories(FTEC_branchlength PUBLIC ${CMAKE_SOURCE_DIR})
