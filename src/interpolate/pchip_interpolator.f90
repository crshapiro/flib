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
module pchip_interpolator
!*******************************************************************************
! The pchip_t class performs linear interpolation for a 1D function v(x). The
! evaluated data is passed on construction. Interpolation is evaluated for real
! sample points xq passed as values or arrays. The first derivatives can also be
! returned as optional arguments.
!
! The resulting interpolation has continuous values and first derivatives.
! Within an interval, the interpolation is monotonic. Extrapolation is used by
! default for evaluation points xq outside of the interval.
!
! For more details see: F. N. Fritsch and R. E. Carlson. (1980). "Monotone
! piecewise cubic interpolation." SIAM Journal on Numerical Analysis. 17(2)
! 238–-246.
!
use stl
use interp1D
implicit none

private
public :: pchip_interpolator_t, pchip_interpolate

type, extends(interp1D_t) :: pchip_interpolator_t
    real(rprec), dimension(:), allocatable :: vp
contains
    procedure, public :: init_pchip_interpolator_t
    procedure, private :: interp_scalar
end type pchip_interpolator_t

interface pchip_interpolator_t
    module procedure constructor
end interface pchip_interpolator_t

interface pchip_interpolate
    module procedure :: pchip_interpolate_scalar
    module procedure :: pchip_interpolate_array
end interface pchip_interpolate

contains

!*******************************************************************************
function pchip_interpolate_scalar(x, v, xq) result(vq)
!*******************************************************************************
! Convenience function interface for pchip_interpolate_t for a single scalar
! query point xq.
!
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
type(pchip_interpolator_t) :: sp

! Create object
sp = pchip_interpolator_t(x,v)

! Interpolate
vq = 0._rprec
call sp%interpolate(xq, vq)

end function pchip_interpolate_scalar

!*******************************************************************************
function pchip_interpolate_array(x, v, xq) result(vq)
!*******************************************************************************
! Convenience function interface for pchip_interpolate_t for an array of
! query points xq.
!
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), allocatable :: vq
type(pchip_interpolator_t) :: sp

! Create object
sp = pchip_interpolator_t(x,v)

! Interpolate
allocate(vq(size(x)))
call sp%interpolate(xq, vq)

end function pchip_interpolate_array

!*******************************************************************************
function constructor(x, v) result(this)
!*******************************************************************************
! Constructor that calls initializer
!
type(pchip_interpolator_t) :: this
real(rprec), dimension(:), intent(in) :: x, v

call this%init_pchip_interpolator_t(x, v)

end function constructor

!*******************************************************************************
subroutine init_pchip_interpolator_t(this, x, v)
!*******************************************************************************
! Initializer for linear_interpolator_t. Takes points v(x) that are used
! for the interpolation. This function also evaluates the second derivative as
! these x's using the tridiagonal matrix algorithm.
!
! By default, the constructor using natural boundary conditions with vanishing
! second derivatives. Other boundary conditions can be specified by supplying
! the optional arguments low_bc or high_bc. The values of the specified
! derivatives can be passed using the optional arguments low_f and high_f.
!
class(pchip_interpolator_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), allocatable :: delta
logical, dimension(:), allocatable :: zeroed
real(rprec) :: a, b, tau
integer :: i

! Call base initializer
call this%init_interp1D_t(x, v, 'extrapolate', 'extrapolate')

! Allocate and set vpp
if ( allocated(this%vp) ) deallocate(this%vp)
allocate( this%vp(this%N) )
this%vp = 0._rprec

! Calculate secant lines between points
allocate( Delta(this%N - 1) )
do i  = 1, this%N-1
    Delta(i) = (this%v(i+1) - this%v(i)) / (this%x(i+1) - this%x(i))
end do

! Set initial guess for the derivatives
!   end points use one-sided differences
!   if the secant lines are different signs, set to zero
!   interior points just average secant lines on both sides
this%vp(1) = Delta(1)
this%vp(this%N) = Delta(this%N-1)
do i = 2, this%N-1
    if ( (Delta(i-1) * Delta(i)) <= 0._rprec ) then
        this%vp(i) = 0._rprec
    else
        this%vp(i) = 0.5_rprec * ( Delta(i-1) + Delta(i) )
    end if
end do

! An array to keep track of zeroed derivative segments
allocate( zeroed(this%N-1) )
zeroed = .false.

! Check if a segment has zero slope and set the derivatives at the end to zero
do i = 1, this%N-1
    if (Delta(i) == 0._rprec) then
        zeroed(i) = .true.
        this%vp(i) = 0._rprec
        this%vp(i+1) = 0._rprec
    end if
end do

! Check that the derivates at the end of each segment have the same sign as the
! secant line. Otherwise, set the derivatives to zero
do i = 1, this%N-1
    if ( .not.zeroed(i) ) then
        a = this%vp(i) / Delta(i)
        b = this%vp(i+1) / Delta(i)
        if ( a < 0._rprec .or. b < 0._rprec ) this%vp(i) = 0._rprec
    end if
end do

! Ensure monotonicity by keeping a and b < 3
do i = 1, this%N-1
    if ( .not.zeroed(i) ) then
        a = this%vp(i) / Delta(i)
        b = this%vp(i+1) / Delta(i)
        tau = 3._rprec / sqrt(a**2 + b**2)
        if (tau < 1._rprec) then
            this%vp(i) = tau * a * Delta(i)
            this%vp(i+1) = tau * b * Delta(i)
        end if
    end if
end do

end subroutine init_pchip_interpolator_t

!*******************************************************************************
subroutine interp_scalar(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for a single point. Uses binary_search to find the
! interval on which the sample point lies. This is a guaranteed log2(N) search
! method.
!
class(pchip_interpolator_t) :: this
real(rprec), intent(in) :: xq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqp
integer :: i
real(rprec) :: t, h
real(rprec) :: a, b, c, d

i = binary_search(this%x, xq)
if (i == 0) then
    vq = this%v(1) + this%vp(1) * (xq - this%x(1))
    if ( present(vqp) ) vqp = this%vp(1)
else if (i == this%N) then
    vq = this%v(this%N) + this%vp(this%N) * (xq - this%x(this%N))
    if ( present(vqp) ) vqp = this%vp(this%N)
else
    h = this%x(i+1)-this%x(i)
    t = ( xq - this%x(i) ) / h
    A = (1._rprec + 2._rprec*t) * (1._rprec - t)**2
    B = t * (1._rprec - t)**2
    C = t**2 * (3._rprec - 2._rprec*t)
    D = t**2 * (t - 1._rprec)
    vq = A*this%v(i) + B*h*this%vp(i) + C*this%v(i+1) + D*h*this%vp(i+1)
    if ( present(vqp) ) then
        A = 6._rprec * t * (t - 1._rprec)
        B = 1._rprec + t * (3._rprec*t - 4._rprec)
        C = 6._rprec * t * (1._rprec - t)
        D = t * (3._rprec*t - 2._rprec)
        vqp = A*this%v(i)/h + B*this%vp(i) + C*this%v(i+1)/h + D*this%vp(i+1)
    end if
end if

end subroutine interp_scalar

end module pchip_interpolator
