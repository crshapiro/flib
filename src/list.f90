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
    type(node_t), pointer, private :: next => null()
    type(node_t), pointer, private :: prev => null()
    class(*), allocatable, private :: item
contains
    final :: node_destructor
end type node_t

interface node_t
    module procedure :: node_constructor
end interface node_t

type :: list_t
    type(node_t), pointer :: front => null()
    type(node_t), pointer :: back => null()
contains
    final :: list_destructor
    procedure, public :: push_front
    procedure, public :: pop_front
    procedure, public :: push_back
    procedure, public :: pop_back
    procedure, public :: insert
    procedure, public :: remove
    procedure, public :: size
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

next => this%front

do while (associated(next))
    this%front => next
    next => this%front%next
    deallocate(this%front)
end do

nullify(this%front)
nullify(this%back)

end subroutine list_destructor

!*******************************************************************************
subroutine push_front(this, item)
!*******************************************************************************
implicit none
class(list_t) :: this
class(*), intent(in) :: item
type(node_t), pointer :: front => null()

if (associated(this%front)) then
    front => this%front
    nullify(this%front)
    allocate(this%front, source=node_t(item, next=front))
    front%prev => this%front
    nullify(front)
else
    allocate(this%front, source=node_t(item))
    this%back => this%front
end if

end subroutine push_front

!*******************************************************************************
function pop_front(this) result(item)
!*******************************************************************************
implicit none
class(list_t) :: this
type(node_t), pointer :: front => null()
class(*), allocatable :: item

if (associated(this%front)) then
    allocate(item, source=this%front%item)
    front => this%front
    this%front => this%front%next
    if (associated(this%front)) then
        this%front%prev => null()
    else
        this%back => null()
    end if
    deallocate(front)
else
    write(*,*) "ERROR: list_t%pop_front called on empty list"
    stop 9
end if

end function pop_front

!*******************************************************************************
subroutine push_back(this, item)
!*******************************************************************************
implicit none
class(list_t) :: this
class(*), intent(in) :: item
type(node_t), pointer :: back => null()

if (associated(this%back)) then
    back => this%back
    nullify(this%back)
    allocate(this%back, source=node_t(item, prev=back))
    back%next => this%back
    nullify(back)
else
    allocate(this%back, source=node_t(item))
    this%front => this%back
end if

end subroutine push_back

!*******************************************************************************
function pop_back(this) result(item)
!*******************************************************************************
implicit none
class(list_t) :: this
type(node_t), pointer :: back => null()
class(*), allocatable :: item

if (associated(this%back)) then
    allocate(item, source=this%back%item)
    back => this%back
    this%back => this%back%prev
    if (associated(this%back)) then
        this%back%next => null()
    else
        this%front => null()
    end if
    deallocate(back)
else
    write(*,*) "ERROR: list_t%pop_back called on empty list"
    stop 9
end if

end function pop_back

!*******************************************************************************
subroutine insert(this, item, n)
!*******************************************************************************
implicit none
class(list_t) :: this
class(*), intent(in) :: item
integer, intent(in) :: n
type(node_t), pointer :: node => null()
type(node_t), pointer :: push_node
integer :: i

if (n <= 0) then
    write(*,*) "ERROR: list_t%insert index must be > 0"
    stop 9
end if

if (n == 1) then
    call this%push_front(item)
    return
end if

node => this%front

do i = 1, n-1
    if (associated(node)) then
        node => node%next
    else
        write(*,*) "ERROR: list_t%insert index greater than list size"
        stop 9
    end if
end do

if (associated(node)) then
    if (associated(node%next)) then
        allocate(push_node, source=node_t(item, prev=node%prev, next=node))
        node%prev%next => push_node
        node%prev => push_node
    else
        call this%push_back(item)
    end if
else
    write(*,*) "ERROR: list_t%insert index greater than list size"
    stop 9
end if

end subroutine insert

!*******************************************************************************
function remove(this, n) result(item)
!*******************************************************************************
implicit none
class(list_t) :: this
integer, intent(in) :: n
class(*), allocatable :: item
type(node_t), pointer :: node => null()
integer :: i

if (n <= 0) then
    write(*,*) "ERROR: list_t%remove index must be > 0"
    stop 9
end if

if (n == 1) then
#ifdef __INTEL_COMPILER
    allocate(item, source=this%pop_front())
#else
    item = this%pop_front()
#endif
    return
end if

node => this%front

do i = 1, n-1
    if (associated(node)) then
        node => node%next
    else
        write(*,*) "ERROR: list_t%remove index greater than list size"
        stop 9
    end if
end do

if (associated(node)) then
    if (associated(node%next)) then
        allocate(item, source=node%item)
        node%prev%next => node%next
        node%next%prev => node%prev
        deallocate(node)
    else
#ifdef __INTEL_COMPILER
        allocate(item, source=this%pop_back())
#else
        item = this%pop_back()
#endif
    end if
else
    write(*,*) "ERROR: list_t%remove index greater than list size"
    stop 9
end if

end function remove

!*******************************************************************************
function size(this) result(n)
!*******************************************************************************
implicit none
class(list_t) :: this
type(node_t), pointer :: node => null()
integer :: n

node => this%front
n = 0

do while (associated(node))
    n = n + 1
    node => node%next
end do

end function size

end module list