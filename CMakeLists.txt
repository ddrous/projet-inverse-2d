cmake_minimum_required(VERSION 3.5)

# Set the project name
project(Inverse)

# Tell CMake to use CTest extension
enable_testing()

# Include the muParser subdirectory
include_directories(src/include)
include_directories(thirdparty/muparser/include)

# Add source files needed
add_executable(transfer 
    src/main.cpp 
    src/exporter.cpp 
    src/solver.cpp 
    src/mesh.cpp 
    src/config.cpp
    thirdparty/muparser/src/muParser.cpp
    thirdparty/muparser/src/muParserBase.cpp
    thirdparty/muparser/src/muParserBytecode.cpp
    thirdparty/muparser/src/muParserCallback.cpp
    thirdparty/muparser/src/muParserError.cpp
    thirdparty/muparser/src/muParserTokenReader.cpp)

# Lots of warnings and all warnings as errors
add_compile_options(-Wall -Wextra -pedantic -Werror)

# A simple test that works perfectly
add_test(NAME simple_run COMMAND transfer ../src/config/test_1.cfg)
# Tests that should fail to make sure the program doesnt crash
add_test(NAME no_config COMMAND transfer)
add_test(NAME bad_param COMMAND transfer ../src/config/test_2.cfg)
add_test(NAME unparsable COMMAND transfer ../src/config/test_3.cfg)
set_tests_properties(no_config bad_param unparsable PROPERTIES WILL_FAIL TRUE)
