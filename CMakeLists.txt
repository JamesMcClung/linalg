cmake_minimum_required(VERSION 3.20.0-rc1)

project(LINALG)

# Set variables

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    set(LINALG_STANDALONE TRUE)
endif()

# Options

option(LINALG_BUILD_TESTS "Build test programs" ${LINALG_STANDALONE})

# Main library

add_library(linalg INTERFACE)
target_include_directories(linalg INTERFACE include/)

# Testing

if (${LINALG_BUILD_TESTS})
    add_subdirectory(tests)
endif()