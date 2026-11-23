cmake_minimum_required(VERSION 3.20)
project(AtomWeave CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(AtomWeave ${CMAKE_SOURCE_DIR}/src/main.cpp)

target_include_directories(AtomWeave
        PRIVATE ${CMAKE_SOURCE_DIR}/src
        PRIVATE ${CMAKE_SOURCE_DIR}/include/atomweave/
)

add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Werror=return-type
)

