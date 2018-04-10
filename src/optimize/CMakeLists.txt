link_libraries(lapack blas)
add_library(optimize lbfgsb.f90 conjugate_gradient.f90 line_search.f90
    minimize.f90 optimize.f90)
target_link_libraries(optimize core)