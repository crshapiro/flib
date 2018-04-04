# Optimization library

## Minimization type (minimize.f90)
**mininize_t** is a base class for use with univariate of multivariate minimization routines. This class can be used as a wrapper for procedures or a base class for other types
### Function wrapper
**minimize_t** can be used to wrap functions with the following abstract interface:
```
subroutine minimize_function(x, f, g)
    use stl
    implicit none
    real(rprec), dimension(:), intent(in) :: x      ! Point to evaluate
    real(rprec), intent(inout) :: f                 ! Function value (scalar)
    real(rprec), dimension(:), intent(inout) :: g   ! Function gradient
end subroutine minimize_function
```
An instance of **minimize_t** wrapping the function *fun* is constructed using the defined constructor:
```
use minimize
implicit none
type(minimize_t) :: m
...
m = minimize_t(fun)
```
### Base class
The **minimize_t** can be used as a base class by overriding the class method *eval*:
```
type, extends(minimize_t) :: derived_t
    ...
contains
    ...
    procedure, public :: eval
    ...
end type derived_t
```
*eval* must have a matching interface:
```
!*******************************************************************************
subroutine eval(this, x, f, g)
!*******************************************************************************
! Evaluates the function pointer. Overload this procedure if extending this base
! class
implicit none
class(derived_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x      ! Point to evaluate
real(rprec), intent(inout) :: f                 ! Function value (scalar)
real(rprec), dimension(:), intent(inout) :: g   ! Function gradient

...

end subroutine eval
```

## Line Search (line_search.f90)

## L-BFGS-B (lbfgsb.f90)

## Conjugate gradient (conjugate_gradient.f90)