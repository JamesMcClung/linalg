# used for testing, as in https://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

enable_testing()

# Order tests by dependency, if possible
add_executable(test_linalg
        test_fullmatrix.cc
        test_zeromatrix.cc
        test_bandedmatrix.cc
        test_recursivematrix.cc
        test_polymorphism.cc
        test_efficiency.cc
        test_io.cc
        test_lu.cc
        test_thomas.cc
        
        test_complex.cc
        test_fft.cc
        )

target_link_libraries(test_linalg linalg gtest_main)

set(TEST_OUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/out/")
configure_file(config.hh.in config.hh)
file(MAKE_DIRECTORY ${TEST_OUT_DIRECTORY})
target_include_directories(test_linalg PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

gtest_discover_tests(test_linalg)

set_target_properties(test_linalg
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/"
)