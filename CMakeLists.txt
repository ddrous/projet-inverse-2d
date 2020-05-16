cmake_minimum_required(VERSION 3.10)

# Go to [https://ncona.com/2019/03/building-a-cpp-project-with-cmake/] for details 
# also [https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html]

# set the project name
project(Inverse)

# tell CMake to use CTest extension
enable_testing()

# Include the exprtk subdirectory
include_directories(thirdparty/exprtk)

# Add source files needed
add_executable(transfer src/main.cpp src/solver.cpp src/mesh.cpp src/config.cpp)

# add_test(spreading spreading) - register the executable with CTest
add_test(NAME simple_run COMMAND transfer ../src/config/test.cfg)
# more tests to come
