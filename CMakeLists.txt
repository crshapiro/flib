cmake_minimum_required(VERSION 2.8.9)
project(flib)
enable_language(Fortran)

enable_testing()

set(CMAKE_Fortran_FLAGS "-cpp -O3 -Wall")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
add_subdirectory(src)