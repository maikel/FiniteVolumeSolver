cmake_minimum_required(VERSION 2.8.12)
project(PackageTest C CXX Fortran)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set (CMAKE_CXX_STANDARD 11)
find_package(AMReX REQUIRED)

add_executable(example example.cpp)
target_link_libraries(example AMReX::amrex)

# CTest is a testing tool that can be used to test your project.
# enable_testing()
# add_test(NAME example
#          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
#          COMMAND example)
