cmake_minimum_required(VERSION 3.0)

project(internal_branches)

add_compile_options(-O3)
#add_compile_options(-ggdb)

include_directories(eigen)
include_directories(tclap)
include_directories(ib_boost)

add_executable(internal_branches internal_branches.cpp)
add_executable(FTEC_branchlength FTEC_branchlength.cpp)
