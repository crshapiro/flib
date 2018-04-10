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
module interpolate
!*******************************************************************************
! The interpolation library currently provides univariate interpolation. Each
! interpolation method can be performed using the classes, where the data v(x)
! are passed on construction. Interpolation is evaluated for real sample points
! xq passed as values or arrays and returned in the variable vq using the
! subroutine *interpolate*.
!
use stl
use linear_interpolator
use cubic_spline_interpolator
use pchip_interpolator
implicit none

contains

!*******************************************************************************
subroutine interpolate_unit_test
!*******************************************************************************
! Interpolation unit test program
!
type(linear_interpolator_t) :: l
type(cubic_spline_interpolator_t) :: cspl
type(pchip_interpolator_t) :: pc
real(rprec), dimension(:), allocatable :: x, v
real(rprec), dimension(:), allocatable :: xq, vq_cspl, vq_pc, vq_l
real(rprec) :: xmin = 0._rprec, xmax = 2._rprec
real(rprec) :: xqmin = -2._rprec, xqmax = 4._rprec
integer :: N = 10, Nq = 100
integer :: i

write(*,*) "*******************************************************************"
write(*,*) " Interpolate unit test"
write(*,*) "*******************************************************************"

! Data points
allocate(x(N), v(N))
do i = 1, N
    x(i) = (i-1)*(xmax-xmin)/(N-1) + xmin
    v(i) = sin(x(i))
end do

! Generate interpolators
l = linear_interpolator_t(x, v)
cspl = cubic_spline_interpolator_t(x, v)
pc = pchip_interpolator_t(x, v)

! Query points
allocate(xq(Nq), vq_cspl(Nq), vq_pc(Nq), vq_l(Nq))
do i = 1, Nq
    xq(i) = (i-1)*(xqmax-xqmin)/(Nq-1) + xqmin
end do
call l%interpolate(xq, vq_l)
call cspl%interpolate(xq, vq_cspl)
call pc%interpolate(xq, vq_pc)

! Output
write(*,*) "                   x                    sin(x)                  "&
    //"  linear              cubic_spline                     pchip"
do i = 1, Nq
    write(*,*) xq(i), sin(xq(i)), vq_l(i), vq_cspl(i), vq_pc(i)
end do

end subroutine interpolate_unit_test

end module interpolate


