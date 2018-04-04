!   Copyright (C) 2017-2018 Carl Shapiro
!
!   This file is part of flib.
!
!   flib is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   flib is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with flib.  If not, see <http://www.gnu.org/licenses/>.

!*******************************************************************************
program optimize
!*******************************************************************************
! Optimization unit test program
use stl
use functions
use minimize
use lbfgsb
use conjugate_gradient
implicit none

type(minimize_t), target :: m
type(lbfgsb_t) :: lb
type(conjugate_gradient_t) :: cg
real(rprec), dimension(:), allocatable :: x

! Initial guess
allocate( x(2) )
x(1) = 50._rprec
x(2) = 12._rprec

! Write information to output
write(*,*) "Demonstrating minimization features using L-BFGS-B with inexact "  &
    // "More-Theunte line search."
write(*,*) "Minimizing the Rosenbrock function..."
write(*,*) "Starting position is (x, y) = ", x

! Perform minimization
m = minimize_t(rosenbrock_wrapper)
lb = lbfgsb_t(m)
call lb%minimize(x)
lb = lbfgsb_t(m)

! Write output
write(*,*) "Optimal point is (x,y) = ", x

! Initial guess
x(1) = 50._rprec
x(2) = 12._rprec

! Write information to output
write(*,*) "Demonstrating minimization features using conjugate gradient with "&
    // "inexact More-Theunte line search."
write(*,*) "Minimizing the Rosenbrock function..."
write(*,*) "Starting position is (x, y) = ", x

! Perform minimization
m = minimize_t(rosenbrock_wrapper)
cg = conjugate_gradient_t(m)
call cg%minimize(x)

! Write output
write(*,*) "Optimal point is (x,y) = ", x

contains

!*******************************************************************************
subroutine rosenbrock_wrapper(x, f, g)
!*******************************************************************************
! Wrapper for Rosenbrock function that conforms to the minimize_function
! abstract interface
implicit none

real(rprec), dimension(:), intent(in) :: x      ! Point to evaluate
real(rprec), intent(inout) :: f                 ! Function value (scalar)
real(rprec), dimension(:), intent(inout) :: g   ! Function gradient

if (size(x) .ne. 2 .or. size(g) .ne. 2) then
    write(*,*) "ERROR: rosenbrock_wrapper: invalid input size"
    stop 9
end if

f = rosenbrock(x(1), x(2))
g(1) = -2*(1._rprec-x(1)) - 400._rprec*(x(2)-x(1)**2)*x(1)
g(2) = 200._rprec*(x(2)-x(1)**2)

end subroutine

end program optimize

