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
module interp1D
!*******************************************************************************
! Abstract base class for 1D interpolators.
!
use stl
implicit none

private
public :: interp1D_t

type, abstract :: interp1D_t
    real(rprec), dimension(:), allocatable :: x, v
    integer :: N
    character(:), allocatable :: low_bc, high_bc
contains
    procedure, public :: init_interp1D_t
    procedure(interp_scalar), private, deferred :: interp_scalar
    procedure, private :: interp_array
    generic, public :: interpolate => interp_scalar, interp_array
end type interp1D_t

abstract interface
    subroutine interp_scalar(this, xq, vq, vqp)
        use stl
        import interp1D_t
        class(interp1D_t) :: this
        real(rprec), intent(in) :: xq
        real(rprec), intent(out) :: vq
        real(rprec), intent(out), optional :: vqp
    end subroutine
end interface

contains

!*******************************************************************************
subroutine init_interp1D_t(this, x, v, low_bc, high_bc)
!*******************************************************************************
! Constructor for interp1D_t. Takes points v(x) that are used for the
! interpolation.
!
class(interp1D_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x, v
character(*), intent(in) :: low_bc, high_bc
integer :: i

! Set the size of the matrix
this%N = size(x)

! Check that all input arguments are the same size
if ( size(v) /= this%N ) then
    write(*,*) "ERROR: interp1D_t%init_interp1D_t: x and v must be the same size"
    stop 9
end if

! Check that x is sorted
do i = 2, this%N
    if ( x(i) < x(i-1) ) then
        write(*,*) "ERROR: interp1D_t%init_interp1D_t: x must be increasing"
        stop 9
    end if
end do

! Deallocate if necessary
if ( allocated(this%x) ) deallocate(this%x)
if ( allocated(this%v) ) deallocate(this%v)

! Allocate and assign points
allocate( this%x(this%N) )
allocate( this%v(this%N) )
this%x = x
this%v = v
this%high_bc = high_bc
this%low_bc = low_bc

end subroutine init_interp1D_t

!*******************************************************************************
subroutine interp_array(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for an array of points. This simply calls interp_scalar
! for each of the sample points.
!
class(interp1D_t) :: this
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), intent(out) :: vq
real(rprec), dimension(:), intent(out), optional :: vqp
integer :: i, N

N = size(xq)

if ( present(vqp) ) then
    do i = 1, N
        call this%interpolate(xq(i), vq(i), vqp(i))
    end do
else
    do i = 1, N
        call this%interpolate(xq(i), vq(i))
    end do
end if

end subroutine interp_array

end module interp1D
