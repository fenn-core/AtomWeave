cmake_minimum_required(VERSION 3.25)

project(AtomWeave LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(atomweave_core INTERFACE)

target_compile_features(atomweave_core INTERFACE cxx_std_20)

target_include_directories(atomweave_core INTERFACE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include/atomweave
)

target_compile_options(atomweave_core INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Werror=return-type
)


add_executable(AtomWeave src/main.cpp)
target_link_libraries(AtomWeave PRIVATE atomweave_core)


include(CTest)

if(BUILD_TESTING)
    add_custom_target(atomweave_tests)

    function(atomweave_add_test name)
        set(target test_${name})

        add_executable(${target} tests/${target}.cpp)

        target_link_libraries(${target} PRIVATE atomweave_core)

        target_include_directories(${target} PRIVATE
                ${PROJECT_SOURCE_DIR}/tests
        )

        add_test(NAME ${name} COMMAND ${target})
        add_dependencies(atomweave_tests ${target})
    endfunction()

    atomweave_add_test(forces)
    atomweave_add_test(dynamics)
endif()
