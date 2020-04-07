cmake_minimum_required(VERSION 3.10)

# set the project name
project(Propagation)

include_directories(src)

add_executable(spreading spreading.cpp solveur.cpp maille.cpp)
add_test(spreading spreading)
