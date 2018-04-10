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
module functions
!*******************************************************************************
! Function library
!
! The function library **functions** contains useful mathematical functions used
! frequently in many contexts
!
use stl
implicit none
private

public :: softplus
interface softplus
    module procedure :: softplus_scalar
    module procedure :: softplus_array
end interface softplus

public :: logistic
interface logistic
    module procedure :: logistic_scalar
    module procedure :: logistic_array
end interface logistic

public :: gaussian
interface gaussian
    module procedure :: gaussian_scalar
    module procedure :: gaussian_array
end interface gaussian

public :: rosenbrock
interface rosenbrock
    module procedure :: rosenbrock_scalar
    module procedure :: rosenbrock_array
end interface rosenbrock

contains

!*******************************************************************************
function softplus_scalar(s, x) result(sp)
!*******************************************************************************
! softplus function of the form sp(x) = ln(1 + exp(x-s)) with scalar input
!
real(rprec), intent(in) :: x, s
real(rprec) :: sp
real(rprec), parameter :: threshold = 100

if (x - s > threshold) then
    sp = x - s
else
    sp = log(1 + exp(x - s))
endif

end function softplus_scalar

!*******************************************************************************
function softplus_array(s, x) result(sp)
!*******************************************************************************
! softplus function of the form ln(1 + exp(x-s)) with array input
!
real(rprec), dimension(:), intent(in)   :: x
real(rprec), intent(in) :: s
real(rprec), dimension(:), allocatable  :: sp
integer :: i

allocate( sp(size(x)) )
do i = 1, size(x)
    sp(i) = softplus(s, x(i))
end do

end function softplus_array

!*******************************************************************************
function logistic_scalar (s, x) result(l)
!*******************************************************************************
! logistic function of the form l(x) = 1/(1 + exp(-(x-s)) with scalar input
!
real(rprec), intent(in) :: x, s
real(rprec) :: l
real(rprec), parameter :: threshold = 100

if (x - s > threshold) then
    l = 1.0
else
    l = 1.0 / ( 1.0 + exp( -(x - s) ) )
endif

end function logistic_scalar

!*******************************************************************************
function logistic_array(s, x) result(l)
!*******************************************************************************
! logistic function of the form l(x) = 1/(1 + exp(-(x-s)) with array input
!
real(rprec), dimension(:), intent(in) :: x
real(rprec), intent(in) :: s
real(rprec), dimension(:), allocatable :: l
integer :: i

allocate( l(size(x)) )
do i = 1, size(x)
    l(i) = logistic(s, x(i))
end do

end function logistic_array

!*******************************************************************************
function gaussian_scalar(x, x0, Delta) result(g)
!*******************************************************************************
! normalized Gaussian with scalar input
!
real(rprec), intent(in) :: x, x0, Delta
real(rprec) :: g

g = 1.0 / (Delta * sqrt(2.0 *pi)) * exp(-0.5 * (x - x0) * (x - x0) / Delta / Delta);
! If near precision limit, set to zero
! NOTE: remove this after check
if (abs(g) < sqrt(epsilon(g))) g = 0.0

end function gaussian_scalar

!*******************************************************************************
function gaussian_array(x, x0, Delta) result(g)
!*******************************************************************************
! normalized Gaussian with array input
!
real(rprec), dimension(:), intent(in) :: x
real(rprec), intent(in) :: x0, Delta
real(rprec), dimension(:), allocatable :: g
integer :: i

allocate( g(size(x)) )
do i = 1, size(x)
    g(i) = gaussian(x(i), x0, Delta)
end do

end function gaussian_array

!*******************************************************************************
function rosenbrock_scalar(x, y, i_a, i_b) result(f)
!*******************************************************************************
! Rosenbrock function f(x,y) = (a-x)^2 + b*(y-x^2)^2 with scalar input
! Default values are a=1 and b=100
!
real(rprec), intent(in) :: x, y
real(rprec), intent(in), optional :: i_a, i_b
real(rprec) :: f
real(rprec) :: a, b

if (present(i_a)) then
    a = i_a
else
    a = 1._rprec
end if

if (present(i_b)) then
    b = i_b
else
    b = 100._rprec
end if

f = (a-x)**2 + b*(y-x**2)**2

end function rosenbrock_scalar

!*******************************************************************************
function rosenbrock_array(x, y, a, b) result(f)
!*******************************************************************************
! Rosenbrock function f(x,y) = (a-x)^2 + b*(y-x^2)^2 with array input
!
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), intent(in), optional :: a, b
real(rprec), dimension(:), allocatable :: f
integer :: i

! Check size and allocate
if (size(x) .ne. size(y)) then
    write(*,*) "ERROR: rosenbrock: x and y must be the same size"
    stop 9
end if
allocate( f(size(x)) )

! Call scalar function, depending on how many optional arguments are included
if (present(a)) then
    if (present(b)) then
        do i = 1, size(x)
            f(i) = rosenbrock(x(i), y(i), a, b)
        end do
    else
        do i = 1, size(x)
            f(i) = rosenbrock(x(i), y(i), a)
        end do
    end if
else
    do i = 1, size(x)
        f(i) = rosenbrock(x(i), y(i))
    end do
end if

end function rosenbrock_array

end module functions
