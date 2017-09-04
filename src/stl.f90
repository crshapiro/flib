!   Copyright (C) 2017 Carl Shapiro
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
function integer_s(cstar) result(i)
!*******************************************************************************
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

end module stl