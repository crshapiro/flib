add_subdirectory(optimize)
add_subdirectory(interpolate)

add_library(core stl.f90 functions.f90 list.f90 tridiagonal.f90)
add_executable(main main.f90)
target_link_libraries(main interpolate optimize)
add_test(main main)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
