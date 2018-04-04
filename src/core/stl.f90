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
module stl
!*******************************************************************************
implicit none

private
public :: uppercase, lowercase, binary_search, count_lines

! rprec is used to specify precision
integer, parameter, public :: rprec = kind(1.d0)
! pi 3.14159.....
real(rprec), parameter, public :: pi = 4._rprec*datan(1._rprec)

interface linear_interp
    module procedure :: linear_interp_ss
    module procedure :: linear_interp_sa
    module procedure :: linear_interp_aa
end interface linear_interp

interface bilinear_interp
    module procedure :: bilinear_interp_ss
    module procedure :: bilinear_interp_sa
    module procedure :: bilinear_interp_aa
end interface bilinear_interp

public :: integer_
interface integer_
    module procedure :: integer_a
    module procedure :: integer_s
end interface integer_

public :: real_
interface real_
    module procedure :: real_s
    module procedure :: real_a
end interface real_

public :: character_
interface character_
    module procedure :: character_a
end interface character_

public :: logical_
interface logical_
    module procedure :: logical_s
    module procedure :: logical_a
end interface logical_

public :: complex_
interface complex_
    module procedure :: complex_s
    module procedure :: complex_a
end interface complex_

contains

!*******************************************************************************
real(rprec) function linear_interp_ss(u1,u2,dx,xdiff)
!*******************************************************************************
!
!  This function performs linear interpolation
!
!  Inputs:
!  u1           - lower bound value in the increasing index direction
!  u2           - upper bound value in the increasing index direction
!  dx           - length delta for the grid in the correct direction
!  xdiff        - distance from the point of interest to the u1 node
!
implicit none

real(rprec), intent(in) :: u1, u2, dx, xdiff

linear_interp_ss = u1 + (xdiff) * (u2 - u1) / dx

end function linear_interp_ss

!*******************************************************************************
function linear_interp_sa_nocheck(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs the linear interpolation for linear_interp_sa
!  and linear_interp_aa without checking bounds of the input arrays.
!  It cannot be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
integer :: i, N

N = size(v)
i = binary_search(x, xq)
if (i == 0) then
    vq = v(1)
else if (i == N) then
    vq = v(N)
else
    vq = linear_interp_ss(v(i), v(i+1), x(i+1)-x(i), (xq - x(i)))
end if

end function linear_interp_sa_nocheck

!*******************************************************************************
function linear_interp_sa(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to a query point xq
!
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq

if ( size(v) /= size(x) ) then
    write(*,*) "ERROR: linear_interp_sa_nocheck: "                             &
    // "Arrays x and v must be the same size"
    stop 9
end if

vq = linear_interp_sa_nocheck(x, v, xq)

end function linear_interp_sa

!*******************************************************************************
function linear_interp_aa(x, v, xq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to an array of query points xq
!
!  Inputs:
!  x  - array of sample points
!  v  - array of values at sample points
!  xq - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N

! Check array sizes
if ( size(v) /= size(x) ) then
    write(*,*) "ERROR: linear_interp_aa: Arrays x and v must be the same size"
    stop 9
end if

! Allocate output
N = size(xq)
allocate(vq(N))

! For each element of the array perform interpolation
do i = 1, N
    vq(i) = linear_interp_sa_nocheck(x, v, xq(i))
end do

end function linear_interp_aa


!*******************************************************************************
real(rprec) function bilinear_interp_ss(u11,u21,u12,u22,dx,dy,xdiff,ydiff)
!*******************************************************************************
!
!  This function performs linear interpolation
!
!  Inputs:
!  u11   - lower bound value in x direction for lower y
!  u21   - upper bound value in x direction for lower y
!  u12   - lower bound value in x direction for upper y
!  u22   - upper bound value in x direction for upper y
!  dx    - length delta for the grid in x direction
!  dy    - length delta for the grid in y direction
!  xdiff - distance from the point of interest to the u11 node in x direction
!  ydiff - distance from the point of interest to the u11 node in y direction
!
implicit none

real(rprec), intent(in) :: u11, u12, u21, u22, dx, dy, xdiff, ydiff
real(rprec) :: v1, v2

v1 = linear_interp(u11, u21, dx, xdiff)
v2 = linear_interp(u12, u22, dx, xdiff)

bilinear_interp_ss = linear_interp(v1,v2,dy,ydiff)

end function bilinear_interp_ss

!*******************************************************************************
function bilinear_interp_sa_nocheck(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs the linear interpolation for bilinear_interp_sa
!  and bilinear_interp_aa without checking bounds of the input arrays.
!  It cannot be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq
integer     :: i, j, Nx, Ny

Nx = size(x)
Ny = size(y)
i = binary_search(x, xq)
if (i == 0) then
    vq = linear_interp(y, v(1,:), yq)
else if (i == Nx) then
    vq = linear_interp(y, v(Nx,:), yq)
else
    j = binary_search(y, yq)
    if (j == 0) then
        vq = linear_interp(v(i,1), v(i+1,1), x(i+1)-x(i), xq - x(i))
    else if (j == Ny) then
        vq = linear_interp(v(i,Ny), v(i+1,Ny), x(i+1)-x(i), xq - x(i))
    else
        vq = bilinear_interp_ss( v(i,j), v(i+1,j), v(i,j+1), v(i+1,j+1), &
             x(i+1) - x(i), y(j+1) - y(j), xq - x(i), yq - y(j) )
    end if
end if

end function bilinear_interp_sa_nocheck

!*******************************************************************************
function bilinear_interp_sa(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points
!  defined on a grid (x,y) with values v to a query point (xq, yq)
!
!  Inputs:
!  x  - array of sample points
!  v  - array of values at sample points
!  xq - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
     write(*,*) "ERROR: bilinear_interp_saa: Array v must have a size of "     &
     // "[size(x), size(y)]"
     stop 9
end if

vq = bilinear_interp_sa_nocheck(x, y, v, xq, yq)

end function bilinear_interp_sa

!*******************************************************************************
function bilinear_interp_aa(x, y, v, xq, yq) result(vq)
!*******************************************************************************
!
!  This function performs linear interpolation from a set of points
!  defined on a grid (x,y) with values v to an array of query points
!  (xq, yq)
!
!  Inputs:
!  x  - array of sample points
!  v  - array of values at sample points
!  xq - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), dimension(:), intent(in) :: xq, yq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
    write(*,*) "ERROR:  bilinear_interp_aa: Array v must have a size of [size(x), size(y)]"
    stop 9
end if
if ( size(xq) /= size(yq) ) then
    write(*,*) "ERROR:  bilinear_interp_aa: Arrays xq and yq must be the same size"
    stop 9
end if

N = size(xq)
allocate(vq(N))

do i = 1, N
    vq(i) = bilinear_interp_sa_nocheck(x, y, v, xq(i), yq(i))
end do

end function bilinear_interp_aa

!*******************************************************************************
function integer_s(cstar) result(i)
!*******************************************************************************
! convert a class(*) scalar to an integer scalar
implicit none
class(*), intent(in) :: cstar
integer :: i

select type (cstar)
    type is (integer)
        i = cstar
    class default
        stop "Invalid argument for integer(class(*))"
end select

end function integer_s

!*******************************************************************************
function integer_a(cstar) result(i)
!*******************************************************************************
! convert a class(*) array to an integer array
implicit none
class(*), dimension(:), intent(in) :: cstar
integer, dimension(:), allocatable :: i

select type (cstar)
    type is (integer)
        allocate(i, source=cstar)
    class default
        stop "Invalid argument for integer(class(*))"
end select

end function integer_a

!*******************************************************************************
function real_s(cstar) result(r)
!*******************************************************************************
! convert a class(*) scalar to a real scalar
implicit none
class(*), intent(in) :: cstar
real :: r

select type (cstar)
    type is (real)
        r = cstar
    class default
        stop "Invalid argument for real(class(*))"
end select

end function real_s

!*******************************************************************************
function real_a(cstar) result(r)
!*******************************************************************************
! convert a class(*) array to a real array
implicit none
class(*), dimension(:), intent(in) :: cstar
real, dimension(:), allocatable :: r

select type (cstar)
    type is (real)
        allocate(r, source=cstar)
    class default
        stop "Invalid argument for real(class(*))"
end select

end function real_a

!*******************************************************************************
function character_a(cstar) result(ch)
!*******************************************************************************
! convert a class(*) array to a character array
implicit none
class(*), intent(in) :: cstar
character(:), allocatable :: ch

select type (cstar)
    type is (character(len=*))
        allocate(ch, source=cstar)
    class default
        stop "Invalid argument for character(class(*))"
end select

end function character_a

!*******************************************************************************
function logical_s(cstar) result(l)
!*******************************************************************************
! convert a class(*) scalar to a logical scalar
implicit none
class(*), intent(in) :: cstar
logical :: l

select type (cstar)
    type is (logical)
        l = cstar
    class default
        stop "Invalid argument for character(class(*))"
end select

end function logical_s

!*******************************************************************************
function logical_a(cstar) result(l)
!*******************************************************************************
! convert a class(*) array to a logical array
implicit none
class(*), dimension(:), intent(in) :: cstar
logical, dimension(:), allocatable :: l

select type (cstar)
    type is (logical)
        allocate(l, source=cstar)
    class default
        stop "Invalid argument for character(class(*))"
end select

end function logical_a

!*******************************************************************************
function complex_s(cstar) result(l)
!*******************************************************************************
! convert a class(*) scalar to a complex scalar
implicit none
class(*), intent(in) :: cstar
complex :: l

select type (cstar)
    type is (complex)
        l = cstar
    class default
        stop "Invalid argument for character(class(*))"
end select

end function complex_s

!*******************************************************************************
function complex_a(cstar) result(l)
!*******************************************************************************
! convert a class(*) array to complex array
implicit none
class(*), dimension(:), intent(in) :: cstar
complex, dimension(:), allocatable :: l

select type (cstar)
    type is (complex)
        allocate(l, source=cstar)
    class default
        stop "Invalid argument for character(class(*))"
end select

end function complex_a


!*******************************************************************************
function uppercase(str) result(ucstr)
!*******************************************************************************
! convert specified string to upper case. Values inside quotation marks are
! ignored
character(*):: str
character(len_trim(str)):: ucstr
integer :: i, ilen, iav, ioffset, iqc, iquote

ilen = len_trim(str)
ioffset = iachar('A')-iachar('a')
iquote = 0
ucstr = str
do i = 1, ilen
    iav = iachar(str(i:i))

    if(iquote==0 .and. (iav==34 .or.iav==39)) then
        iquote = 1
        iqc = iav
        cycle
    end if

    if(iquote==1 .and. iav==iqc) then
        iquote = 0
        cycle
    end if

    if (iquote==1) cycle

    if(iav >= iachar('a') .and. iav <= iachar('z')) then
        ucstr(i:i) = achar(iav+ioffset)
    else
        ucstr(i:i) = str(i:i)
    end if
end do

end function uppercase

!*******************************************************************************
function lowercase(str) result(lcstr)
!*******************************************************************************
! Convert specified string to lower case. Values inside quotation marks are
! ignored
character(*):: str
character(len_trim(str)):: lcstr
integer :: i, ilen, iav, ioffset, iqc, iquote

ilen = len_trim(str)
ioffset = iachar('A')-iachar('a')
iquote = 0
lcstr = str
do i = 1, ilen
    iav = iachar(str(i:i))

    if(iquote==0 .and. (iav==34 .or.iav==39)) then
        iquote = 1
        iqc = iav
        cycle
    end if

    if(iquote==1 .and. iav==iqc) then
        iquote = 0
        cycle
    end if

    if (iquote==1) cycle

    if(iav >= iachar('A') .and. iav <= iachar('Z')) then
        lcstr(i:i) = achar(iav-ioffset)
    else
        lcstr(i:i) = str(i:i)
    end if
end do

end function lowercase

!*******************************************************************************
function binary_search(arr, val) result(low)
!*******************************************************************************
! Perform a binary search on a sorted array. Given the  provided value,
! adjacent low and high indices of the array are found
! such that the provided value is bracketed. Guaranteed log2(N) search.
!
!  Inputs:
!  arr          - sorted array of values to search
!  val          - value to be bracketed
!
!  Output:
!  low          - lower index of the array bracket
!                 0 if val < arr(1), N if val < arr(N))
!
implicit none

real(rprec), dimension(:) :: arr
real(rprec) :: val
integer :: low, mid, high, N

! Size of array
N = size(arr)

! Check if value is outside bounds
if ( val < arr(1) ) then
    low = 0
    return
end if
if ( val > arr(N) ) then
    low = N
    return
end if

! Otherwise perform bisection
low = 1
high = N
do while (high - low > 1)
    mid = (low + high) / 2
    if ( arr(mid) > val ) then
        high = mid
    elseif ( arr(mid) < val ) then
        low = mid
    else
        low = mid
        return
    endif
end do

end function binary_search

!*******************************************************************************
function count_lines(fname) result(N)
!*******************************************************************************
! Counts the number of lines in a file
implicit none
character(*), intent(in) :: fname
logical :: exst
integer :: fid, ios
integer :: N

! Check if file exists and open
inquire (file=fname, exist=exst)
if (exst) then
    open(newunit=fid, file=fname, position='rewind', form='formatted')
else
    write(*,*) "ERROR: count_lines: file " // trim(fname) // " does not exist"
    stop 9
end if

! count number of lines and close
ios = 0
N = 0
do
    read(fid, *, IOstat = ios)
    if (ios /= 0) exit
    N = N + 1
end do

! Close file
close(fid)

end function count_lines

end module stl
