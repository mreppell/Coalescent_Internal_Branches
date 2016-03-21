cmake_minimum_required(VERSION 3.0)

project(internal_branches)

add_compile_options(-O3)
include_directories(eigen)
include_directories(tclap)
include_directories(boost)

add_executable(internal_branches internal_branches.cpp)
add_executable(FTEC_branchlength FTEC_branchlength.cpp)
