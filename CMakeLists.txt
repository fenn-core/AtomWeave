cmake_minimum_required(VERSION 3.20)
project(AtomWeave LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Werror=return-type
)

add_executable(AtomWeave src/main.cpp)

target_include_directories(AtomWeave PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/include/atomweave
)

enable_testing()

add_executable(test_forces tests/test_dynamics.cpp)
target_include_directories(test_forces PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/include/atomweave
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
)

add_executable(test_dynamics tests/test_dynamics.cpp)
target_include_directories(test_dynamics PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/include/atomweave
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
)

add_test(NAME forces COMMAND test_forces)
add_test(NAME dynamics COMMAND test_dynamics)
