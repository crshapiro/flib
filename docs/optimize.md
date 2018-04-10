# Optimization library

## Minimization type
**mininize_t** is a base class for use with univariate or multivariate minimization routines. This class can be used as a wrapper for procedures or a base class for other types

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

## Line Search
**line_search_t** performs an inexact line search for multivariable minimization. The Algorithm described in J.J. More and D.J. Thuente (1994) "Line search algorithms with guaranteed sufficient decrease." ACM Trans. Math. Software, 20(3), 286--307.

The original source code was written in F77 at: http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/

## L-BFGS-B
**lbrgsb_t** minimizes a function using the L-BFGS-B algorithm.

See J. Nocedal and J.L. Morales (2011) "Remark on "Algorithm 778:
L-BFGS-B: Fortran subroutines for large-scale bound constrained
optimization"  (2011). ACM Transactions on Mathematical Software 38(1)

## Conjugate gradient (conjugate_gradient.f90)
**conjugate_gradient_t** minimizes a function using the conjugate gradient algorithm.