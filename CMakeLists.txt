cmake_minimum_required(VERSION 2.8.5)

# Go to [https://ncona.com/2019/03/building-a-cpp-project-with-cmake/] for details 

# set the project name
project(Propagation)

# Add source files needed
add_executable(Propagation src/spreading.cpp src/solveur.cpp src/maille.cpp src/config.cpp)

# Add the Eigen Librady
add_subdirectory(thirdparty/eigen-3.3.7)

# Linking the Eigen and Propagation
target_link_libraries(Propagation Eigen3)

# add_test(spreading spreading)
