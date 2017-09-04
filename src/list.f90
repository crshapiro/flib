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
module list
!*******************************************************************************
use stl
implicit none

private
public list_t

type :: node_t
    type(node_t), pointer :: next => null()
    type(node_t), pointer :: prev => null()
    class(*), allocatable :: item
contains
    final :: node_destructor
    procedure, public :: destroy
end type node_t

interface node_t
    module procedure :: node_constructor
end interface node_t

type :: list_t
    type(node_t), pointer :: front => null()
    type(node_t), pointer :: back => null()
    integer :: N = 0
contains
    final :: list_destructor
    procedure, public :: push_front
    procedure, public :: pop_front
end type list_t

interface list_t
    module procedure :: list_constructor
end interface list_t

contains

!*******************************************************************************
function node_constructor(item, next, prev) result(this)
!*******************************************************************************
implicit none
class(*), intent(in) :: item
type(node_t), intent(in), pointer, optional :: next, prev
type(node_t) :: this

allocate(this%item, source=item)
if (present(next)) this%next => next
if (present(prev)) this%prev => prev

end function node_constructor

!*******************************************************************************
subroutine node_destructor(this)
!*******************************************************************************
implicit none
type(node_t) :: this

if (allocated(this%item)) deallocate(this%item)
if (associated(this%prev)) nullify(this%prev)
if (associated(this%next)) nullify(this%next)

end subroutine node_destructor


!*******************************************************************************
function list_constructor() result(this)
!*******************************************************************************
implicit none
type(list_t) :: this

end function list_constructor

!*******************************************************************************
subroutine list_destructor(this)
!*******************************************************************************
implicit none
type(list_t) :: this
type(node_t), pointer :: next => null()
integer :: i

write(*,*) "IN LIST_DESTRUCTOR"

next => this%front

do i = 1, this%N
    write(*,*) "item", i
    this%front => next
    next => this%front%next
    deallocate(this%front)
end do

nullify(this%front)

end subroutine list_destructor

!*******************************************************************************
recursive subroutine destroy(this)
!*******************************************************************************
class(node_t) :: this
type(node_t), pointer :: next => null()

write(*,*) "IN DESTROY_NODE"

next => this%next

if (allocated(this%item)) deallocate(this%item)
if (associated(this%prev)) nullify(this%prev)
if (associated(this%next)) nullify(this%next)
if (associated(next)) call next%destroy

end subroutine destroy

!*******************************************************************************
subroutine push_front(this, item)
!*******************************************************************************
implicit none
class(list_t) :: this
class(*), intent(in) :: item
type(node_t), allocatable, target :: push_node
type(node_t), pointer :: front => null()

if (this%N == 0) then
    allocate(this%front, source=node_t(item))
    this%back => this%front
else
    front => this%front
    nullify(this%front)
    allocate(this%front, source=node_t(item, next=front))
    front%prev => this%front
    nullify(front)
end if

this%N = this%N + 1

end subroutine push_front

!*******************************************************************************
function pop_front(this) result(item)
!*******************************************************************************
implicit none
class(list_t) :: this
type(node_t), pointer :: front => null()
class(*), allocatable :: item

if (this%N == 0) then
    write(*,*) "ERROR: list_t%pop_front called on empty list"
    stop 9
else if (this%N == 1) then
    allocate(item, source=this%front%item)
    deallocate(this%front)
    this%front => null()
    this%back => null()
else
    allocate(item, source=this%front%item)
    front => this%front
    this%front => this%front%next
    this%front%prev => null()
    deallocate(front)
end if

this%N = this%N - 1

end function pop_front

end module list