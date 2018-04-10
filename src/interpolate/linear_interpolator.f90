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
module linear_interpolator
!*******************************************************************************
! The linear_interpolator_t class performs linear interpolation for a 1D
! function v(x). The evaluated data is passed on construction.
! Interpolation is evaluated for real sample points xq passed as values or
! arrays. The first derivatives can also be returned as optional arguments.
!
! The resulting interpolation has continuous values and discontinuous first
! derivatives. Extrapolation is used by default for evaluation points xq outside
! of the interval, but extrapolation can be disabled on construction.
!
use stl
use interp1D
implicit none

private
public :: linear_interpolator_t, linear_interpolate

type, extends(interp1D_t) :: linear_interpolator_t
    real(rprec), dimension(:), allocatable :: vp
contains
    procedure, public :: init_linear_interpolator_t
    procedure, private :: interp_scalar
end type linear_interpolator_t

interface linear_interpolator_t
    module procedure constructor
end interface linear_interpolator_t

interface linear_interpolate
    module procedure :: linear_interpolate_scalar
    module procedure :: linear_interpolate_array
end interface linear_interpolate

contains

!*******************************************************************************
function linear_interpolate_scalar(x, v, xq, low_bc, high_bc) result(vq)
!*******************************************************************************
! Convenience function interface for linear_interpolator_t for a single
! scalar query point xq.
!
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
character(*), intent(in), optional :: low_bc, high_bc
real(rprec) :: vq
type(linear_interpolator_t) :: l

! Create object
if (present(low_bc)) then
    if (present(high_bc)) then
        l = linear_interpolator_t(x, v, low_bc, high_bc)
    else
        l = linear_interpolator_t(x, v, low_bc)
    end if
else
    l = linear_interpolator_t(x, v)
end if

! Interpolate
vq = 0._rprec
call l%interpolate(xq, vq)

end function linear_interpolate_scalar

!*******************************************************************************
function linear_interpolate_array(x, v, xq, low_bc, high_bc) result(vq)
!*******************************************************************************
! Convenience function interface for linear_interpolator_t for an array of
! query points xq.
!
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), intent(in) :: xq
character(*), intent(in), optional :: low_bc, high_bc
real(rprec), dimension(:), allocatable :: vq
type(linear_interpolator_t) :: l

! Create object
if (present(low_bc)) then
    if (present(high_bc)) then
        l = linear_interpolator_t(x, v, low_bc, high_bc)
    else
        l = linear_interpolator_t(x, v, low_bc)
    end if
else
    l = linear_interpolator_t(x, v)
end if

! Interpolate
allocate(vq(size(x)))
call l%interpolate(xq, vq)

end function linear_interpolate_array

!*******************************************************************************
function constructor(x, v, i_low_bc, i_high_bc) result(this)
!*******************************************************************************
! Constructor that calls initializer
!
type(linear_interpolator_t) :: this
real(rprec), dimension(:), intent(in) :: x, v
character(*), intent(in), optional :: i_low_bc, i_high_bc
character(:), allocatable :: low_bc, high_bc

! Set defaults
low_bc = "extrapolate"
high_bc = "extrapolate"

! Check optional arguments
if ( present(i_low_bc) ) low_bc = i_low_bc
if ( present(i_high_bc) ) high_bc = i_high_bc

! create object
call this%init_linear_interpolator_t(x, v, low_bc, high_bc)

end function constructor

!*******************************************************************************
subroutine init_linear_interpolator_t(this, x, v, low_bc, high_bc)
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
class(linear_interpolator_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x, v
character(*), intent(in) :: low_bc, high_bc
integer :: i

! Call base initializer
call this%init_interp1D_t(x, v, low_bc, high_bc)

! Check boundary values
select case (uppercase(this%low_bc))
    case('EXTRAPOLATE')
    case('BOUNDARY')
    case default
        write(*,*) "ERROR: linear_interpolator_t: " //                     &
            "Invalid low boundary condition type " // this%low_bc
    stop 9
end select

select case (uppercase(this%high_bc))
    case('EXTRAPOLATE')
    case('BOUNDARY')
    case default
        write(*,*) "ERROR: linear_interpolator_t: " //                     &
            "Invalid high boundary condition type " // this%high_bc
    stop 9
end select

! Allocate and set derivatives
if ( allocated(this%vp) ) deallocate(this%vp)
allocate(this%vp(this%N-1))
do i = 1, this%N-1
    this%vp(i) = (this%v(i+1) - this%v(i)) / (this%x(i+1) - this%x(i))
end do

end subroutine init_linear_interpolator_t

!*******************************************************************************
subroutine interp_scalar(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for a single point. Uses binary_search to find the
! interval on which the sample point lies. This is a guaranteed log2(N) search
! method.
!
class(linear_interpolator_t) :: this
real(rprec), intent(in) :: xq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqp
integer :: i

i = binary_search(this%x, xq)
if (i == 0) then
    select case (uppercase(this%low_bc))
        case('EXTRAPOLATE')
            vq = this%v(1) + (xq - this%x(1)) * this%vp(1)
            if (present(vqp)) vqp = this%vp(1)
        case('BOUNDARY')
            vq = this%v(1)
            if (present(vqp)) vqp = 0._rprec
    end select
else if (i == this%N) then
    select case (uppercase(this%high_bc))
        case('EXTRAPOLATE')
            vq = this%v(this%N) + (xq - this%x(this%N)) * this%vp(this%N-1)
            if (present(vqp)) vqp = this%vp(this%N-1)
        case('BOUNDARY')
            vq = this%v(this%N)
            if (present(vqp)) vqp = 0._rprec
    end select
else
    vq = this%v(i) + (xq - this%x(i)) * this%vp(i)
    if (present(vqp)) vqp = this%vp(i)
end if

end subroutine interp_scalar

end module linear_interpolator
