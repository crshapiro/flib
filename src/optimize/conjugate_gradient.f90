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
module conjugate_gradient
!*******************************************************************************
! This module minimizes a function using the Polak-Ribiere conjugate gradient
! algorithm.
!
use stl
use line_search
use minimize
implicit none

private

public :: conjugate_gradient_t
type :: conjugate_gradient_t
    class(minimize_t), pointer :: mini => NULL()
    ! A small number
    real(rprec) :: eps = 1E-10
    ! maximum number of CG iterations
    integer :: maxiter = 10000
    ! convergence level
    real(rprec) :: tol = 1E-6
    ! current function evaluation
    real(rprec) :: f = 0.0
    ! previous function evaluation
    real(rprec) :: fp = 0.0
    ! number of function evaluations
    integer :: fnev = 0
    ! conjugate direction step
    real(rprec) :: gamma = 0.0
    ! line search class
    type(line_search_t) :: ls
    ! conjugate direction
    real(rprec), dimension(:), allocatable :: gd
    ! current and previous location
    real(rprec), dimension(:), allocatable :: x, xp
    ! current and previous gradients
    real(rprec), dimension(:), allocatable :: g, gp
    ! current and previous search direction
    real(rprec), dimension(:), allocatable :: h, hp
    ! array of lower bounds
    real(rprec), dimension(:), allocatable :: lb
    ! array of upper bounds
    real(rprec), dimension(:), allocatable :: ub
contains
   procedure, public :: minimize
   procedure, private :: evaluate_gamma
end type conjugate_gradient_t

interface conjugate_gradient_t
    module procedure :: constructor
end interface conjugate_gradient_t

contains

!*******************************************************************************
function constructor(i_mini, i_maxiter, i_lb, i_ub, i_tol) result(this)
!*******************************************************************************
! Constructor for conjugate gradient that takes as an argument a pointer to a
! minimize_t
!
type(conjugate_gradient_t) :: this
class(minimize_t), target :: i_mini
integer, intent(in), optional :: i_maxiter
real(rprec), intent(in), dimension(:), optional :: i_lb, i_ub
real(rprec), intent(in), optional :: i_tol

! Assign input arguments
if ( present(i_maxiter) ) this%maxiter = i_maxiter
if ( present(i_tol) ) this%tol = i_tol
if ( present(i_lb) ) then
    allocate(this%lb(size(i_lb)))
    this%lb = i_lb
end if
if ( present(i_ub) ) then
    allocate(this%ub(size(i_ub)))
    this%ub = i_ub
end if
this%mini => i_mini
this%ls = line_search_t(i_mini)

end function constructor

!*******************************************************************************
subroutine evaluate_gamma(this)
!*******************************************************************************
! Calculate the gamma value for the search direction. Uses Polak_Ribiere.
!
class(conjugate_gradient_t), intent(inout) :: this

this%gd = this%g - this%gp
this%gamma = sum(this%g * this%gd) / sum(this%gp * this%gp)

end subroutine evaluate_gamma

!*******************************************************************************
subroutine minimize(this, x)
!*******************************************************************************
! Minimize the function and return the result in the array x
!
class(conjugate_gradient_t), intent(inout) :: this
real(rprec), dimension(:), intent(inout) :: x
real(rprec) :: d, delta_f, stp
integer :: i, j
integer :: dummy = 0

! Allocate arrays
allocate(this%xp(size(x)))
allocate(this%gp(size(x)))
allocate(this%hp(size(x)))
allocate(this%x(size(x)))
allocate(this%g(size(x)))
allocate(this%h(size(x)))
allocate(this%gd(size(x)))

! Evaluate gradient and function at starting position
this%xp = x
call this%mini%eval(this%xp, this%fp, this%gp)
write(*,*) this%fp
this%fnev = this%fnev + 1
this%hp = -this%gp

! Assign other arrays
this%f  = this%fp
this%g  = this%gp
this%h  = this%hp
this%x  = this%xp
this%gd = this%g

d = sum(this%g * this%h)
delta_f = -0.5 * d

do i = 1, this%maxiter
    ! Compute derivative with respect to stp at origin
    d = sum(this%g * this%h)

    ! If derivative is positive, reset at the gradient
    if (d > 0) then
        this%h = -this%g
        d = sum(this%g * this%h)
        delta_f = -0.5 * d
    end if

    ! Search in the direction hp
    stp = min(1.0, -2.0 * delta_f / d)
    call this%ls%search(this%x, this%f, this%g, this%h, stp, dummy)
    this%fnev = this%fnev + dummy

    ! Safeguard step against bounds
    do j = 1, size(this%x)
        if (allocated(this%lb)) this%x(j) = max(this%x(j), this%lb(j))
        if (allocated(this%ub)) this%x(j) = min(this%x(j), this%ub(j))
    end do
    call this%mini%eval(this%x, this%f, this%g)
    write(*,*) this%f

    ! Check for convergence
    if ( 2.0 * abs(this%fp - this%f) <= this%tol * ( abs(this%fp) +            &
    abs(this%f) + this%eps ) .or. sum(this%g*this%g) == 0 ) then
        ! Set output if present
        x = this%x

        ! Evaluate minimization at current point
        call this%mini%eval(this%x, this%f, this%g)

        ! Print result
        write(*,*) 'Conjugate gradient terminated after ', i,                  &
            'iterations. Minimum f = ', this%f

        return
    else if (this%f > this%fp) then
        ! Reset to previous point
        this%x = this%xp;
        this%f = this%fp;
        this%g = this%gp;

        ! Set output if present
        x = this%x

        ! Evaluate minimization at current point
        call this%mini%eval(this%x, this%f, this%g)

        ! Print result
        write(*,*) 'Conjugate gradient terminated after ', i,                  &
            'iterations. Minimum f = ', this%f
        return
    end if

    ! Compute next search direction
    call this%evaluate_gamma
    this%h = -this%g + this%gamma * this%hp

    ! Swap arrays
    delta_f = max(this%fp - this%f, 10 * epsilon(this%f))

    this%gp = this%g
    this%hp = this%h
    this%xp = this%x
    this%fp = this%f
end do

! Evaluate minimization at current point
call this%mini%eval(this%x, this%f, this%g)

! Set output if present
x = this%x

! Print result
write(*,*) 'Conjugate gradient terminated after ', this%maxiter,               &
    'iterations. Minimum f = ', this%f

end subroutine minimize

end module conjugate_gradient
