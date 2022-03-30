cmake_minimum_required(VERSION 3.20.0-rc1)

# used to facilitate linting; see https://stackoverflow.com/questions/36122358/visual-studio-code-c-include-path
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

project(LINALG)

add_library(linalg INTERFACE)
target_include_directories(linalg INTERFACE include/)

add_subdirectory(tests)