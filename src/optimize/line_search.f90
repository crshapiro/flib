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
module line_search
!*******************************************************************************
! Inexact line search for multivariable minimization. Algorithm described in
! J.J. More and D.J. Thuente "Line Search Algorithms with Guaranteed Sufficient
! Decrease." ACM Trans. Math. Software, 20(3), 286--307, 1994. Original source
! code written in F77 at: http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/
use stl
use minimize
implicit none

private
public line_search_t

type :: line_search_t
    ! Pointer to a minimize class
    class(minimize_t), pointer :: mini => NULL()
    ! Sufficient decrease condition parameter
    real(rprec) :: ftol = 0.1
    ! Curvature condition parameter
    real(rprec) :: gtol = 0.4
    ! Maximum relative width of uncertainty
    real(rprec) :: xtol = 1E-6
    ! Minimum step size
    real(rprec) :: minStep = 0.0
    ! Maximum step size
    real(rprec) :: maxStep = 10000
    ! Maximum number of function evaluations
    integer :: maxFunEval = 100
    ! The minimum is currently bracketed
    logical :: bracketed = .false.
    ! Minimum step for current trial within interval of uncertainty
    real(rprec) :: curMinStep  = 0.0
    ! Maximum step for current trial within interval of uncertainty
    real(rprec) :: curMaxStep  = 0.0
    ! Output search details
    logical :: verbose
contains
    ! Public method to search
    procedure, public :: search
    ! Takes a single step
    procedure, private :: take_step
end type line_search_t

interface line_search_t
    module procedure :: constructor
end interface line_search_t

contains

!*******************************************************************************
function constructor(i_mini, i_verbose, i_ftol, i_gtol, i_xtol, i_minStep,     &
    i_maxStep, i_maxFunEval) result(this)
!*******************************************************************************
implicit none
type(line_search_t) :: this
class(minimize_t), target :: i_mini
logical, intent(in), optional :: i_verbose
real(rprec), intent(in), optional :: i_ftol, i_gtol, i_xtol, i_minStep
real(rprec), intent(in), optional :: i_maxStep
integer, intent(in), optional :: i_maxFunEval

! Assign input arguments
this%mini => i_mini
if ( present(i_verbose) ) then
    this%verbose = i_verbose
end if
if ( present(i_ftol) ) then
    this%ftol = i_ftol
end if
if ( present(i_gtol) ) then
    this%gtol = i_gtol
end if
if ( present(i_xtol) ) then
    this%xtol = i_xtol
end if
if ( present(i_minStep) ) then
    this%minStep = i_minStep
end if
if ( present(i_maxStep) ) then
    this%maxStep = i_maxStep
end if
if ( present(i_maxFunEval) ) then
    this%maxFunEval = i_maxFunEval
end if

end function constructor

!*******************************************************************************
subroutine search(this, x, f, g, s, stp, o_numFunEval)
!*******************************************************************************
implicit none
class(line_search_t), intent(inout) :: this
real(rprec), dimension(:), intent(inout) :: x
real(rprec), dimension(:), intent(inout) :: g
real(rprec), intent(inout) :: f, stp
real(rprec), dimension(:), intent(in) :: s
integer, intent(out), optional :: o_numFunEval
integer :: numFunEval
real(rprec), dimension(:), allocatable :: x0
logical :: take_step_success = .true.
logical :: stage1 = .true.
real(rprec) :: f0, dd0, ddtest, width, width1
real(rprec) :: stx, fx, ddx
real(rprec) :: sty, fy, ddy
real(rprec) :: dg, ftest1
real(rprec) :: fm, fxm, fym, dgm, ddxm, ddym

! Initialization
numFunEval = 0
this%bracketed = .false.

! Check input arguments
if (size(x) == 0) then
    write(*,*) "ERROR: line_search_t%search: initial guess must be non-empty."
    stop 9
end if
if (size(x) /= size(g)) then
    write(*,*) "ERROR: line_search_t%search: "                                 &
        // "gradient and input must be the same size."
    stop 9
end if
if (size(x) /= size(s)) then
    write(*,*) "ERROR: line_search_t%search: x and s must be the same size."
    stop 9
end if
if (stp <= 0.0) then
    write(*,*) "ERROR: line_search_t%search: "                                 &
        // "initial step size guess must be > 0."
    stop 9
end if
if (this%minStep < 0.0) then
    write(*,*) "ERROR: line_search_t%search: "                                 &
        // "maximum step size must be greater than minimum step size."
    stop 9
end if
if (this%maxStep < this%minStep) then
    write(*,*) "ERROR: line_search_t%search: "                                 &
        // "minimum step size must be non-negative."
    stop 9
end if
if (this%ftol < 0.0 .or. this%gtol < 0.0 .or. this%xtol < 0.0) then
    write(*,*) "ERROR: line_search_t%search: search tolerances must be >= 0."
    stop 9
end if

! Initialization point x0 and associated function evaluation f0 and
! directional derivative dd0
allocate(x0(size(x)))
x0 = x
f0 = f
dd0 = sum(g*s)
if (dd0 >= 0.0) then
    write(*,*) "ERROR: line_search_t%search: gradient is non-decreasing."
    stop 9
end if

! Initialize local variables.
ddtest = this%ftol * dd0
width  = this%maxStep -  this%minStep
width1 = 2*width

! step, function, and directional derivative at the best step.
stx = 0.0
fx = f0
ddx = dd0
! step, function, and derivative at the other endpoint of the interval
! of uncertainty.
sty = 0.0
fy = f0
ddy = dd0
! Start of iteration
do
    ! Interval of uncertainty
    if (this%bracketed) then
        this%curMinStep = min(stx, sty)
        this%curMaxStep = max(stx, sty)
    else
        this%curMinStep = stx
        this%curMaxStep = stp + 4.0 * (stp-stx)
    end if

    ! Force the step to be within the bounds maxStep and minStep.
    stp = max(stp, this%minStep);
    stp = min(stp, this%maxStep);
    ! If an unusual termination is to occur then let stp be the lowest point
    ! obtained so far.
    if ( (this%bracketed .and. (stp <= this%curMinStep                         &
    .or. stp >= this%curMaxStep)) .or. numFunEval >= this%maxFunEval - 1       &
    .or. .not.take_step_success .or.  (this%bracketed .and. this%curMaxStep    &
    - this%curMinStep <= this%xtol * this%curMaxStep) )  then
        stp = stx
    end if

    ! Evaluate the function and gradient at stp and compute the directional
    ! derivative.
    x = x0 + stp * s
    call this%mini%eval(x, f, g)
    numFunEval = numFunEval + 1
    dg = sum(g*s)
    ftest1 = f0 + stp * ddtest

    ! Tests for convergence issues. Strange behaviors will return normally and
    ! print a warning to std output
    if ((this%bracketed .and. (stp <= this%curMinStep                          &
    .or. stp >= this%curMaxStep))  .or. .not.take_step_success) then
        if (this%verbose) write(*,*) 'line_search_t%search: Rounding errors '  &
            // 'prevent further progress. There may not be a step which '      &
            // 'satisfies the sufficient decrease and curvature conditions. '  &
            // 'Tolerances may be too small.'
        o_numFunEval = numFunEval
        return
    end if
    if (stp == this%maxStep .and. f <= ftest1 .and. dg < ddtest) then
        If (this%verbose) write(*,*) 'line_search_t%search: The step is at '   &
        // 'the upper bound.'
        o_numFunEval = numFunEval
        return
    end if
    if (stp == this%minStep .and. (f > ftest1 .or. dg >= ddtest)) then
        if (this%verbose) write(*,*) 'line_search_t%search: The step is at the'&
        // 'lower bound.'
        o_numFunEval = numFunEval
        return
    end if
    if (numFunEval >= this%maxFunEval) then
        if (this%verbose) write(*,*) 'line_search_t%search: Maximum number of '&
        // 'function evaluations has been reached.'
        o_numFunEval = numFunEval
        return
    end if
    if (this%bracketed .and. this%curMaxStep - this%curMinStep <= this%xtol    &
    * this%curMaxStep) then
        if (this%verbose) write(*,*) 'line_search_t%search: Relative width of '&
        // 'the interval of uncertainty has reached tolerance.'
        o_numFunEval = numFunEval
        return
    end if

    ! Test for convergence. i.e. the sufficient decrease and curvature
    ! conditions hold.
    if ( f <= ftest1 .and. abs(dg) <= this%gtol * (-dd0) ) then
       o_numFunEval = numFunEval
       return
    end if

    ! In the first stage we seek a step for which the modified function has a
    ! nonpositive value and nonnegative
    if ( stage1 .and. f <= ftest1 .and. dg >= min(this%ftol, this%gtol)        &
    * dd0 ) then
        stage1 = .false.
    end if

    ! A modified function is used to predict the step only if we have not
    ! obtained astep for which the modified function has a nonpositive function
    ! value andnonnegative derivative, and if a lower function value has been
    ! obtained but thedecrease is not sufficient.
    if (stage1 .and. f <= fx .and. f > ftest1) then
        ! Define the modified function and derivative values.
        fm = f - stp * ddtest
        fxm = fx - stx * ddtest
        fym = fy - sty * ddtest
        dgm = dg - ddtest
        ddxm = ddx - ddtest
        ddym = ddy - ddtest

        ! Call cstep to update the interval of uncertainty and to compute
        ! the new step.
        call this%take_step(stx, fxm, ddxm, sty, fym, ddym, stp, fm, dgm,      &
            take_step_success)

        ! Reset the function and gradient values for f.
        fx = fxm + stx * ddtest
        fy = fym + sty * ddtest
        ddx = ddxm + ddtest
        ddy = ddym + ddtest
    else
        ! Call cstep to update the interval of uncertainty and to compute the
        ! new step.
        call this%take_step(stx, fx, ddx, sty, fy, ddy, stp, f, dg,            &
            take_step_success)
    end if

    ! Force a sufficient decrease in the size of the interval of uncertainty.
    if (this%bracketed) then
        if (abs(sty-stx) >= 2.0/3.0 * width1) then
            stp = stx + 0.5 * (sty - stx)
        end if
        width1 = width
        width = abs(sty - stx)
    end if

    continue
end do

end subroutine search

!*******************************************************************************
subroutine take_step(this, stx, fx, dx, sty, fy, dy, stp, fp, dp,              &
    take_step_success)
!*******************************************************************************
implicit none
class(line_search_t), intent(inout) :: this
real(rprec), intent(inout) :: stx, fx, dx, sty, fy, dy, stp, fp, dp
logical, intent(inout)     :: take_step_success
integer :: info = 0
real(rprec) :: sgnd
logical :: bound
real(rprec) :: theta, s, gamma
real(rprec) :: p, q, r, stpc, stpq, stpf

! Check the input parameters for errors.
if ( (this%bracketed .and. (stp <= min(stx, sty) .or. stp >= max(stx, sty)))   &
 .or. dx * (stp - stx) >= 0.0 .or. this%curMaxStep < this%curMinStep) then
    write(*,*) "ERROR: line_search%take_step: Invalid input argument"
    stop 9
end if

! Determine if the derivatives have opposite sign.
sgnd = dp * ( dx / abs(dx) );

! First case. A higher function value. The minimum is bracketed. If the cubic
! step is closer to stx than the quadratic step, the cubic step is taken, else
! the average of the cubic and quadratic steps is taken.
if (fp > fx) then
    info = 1
    bound = .true.
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
    s = max(abs(theta), max(abs(dx),abs(dp)))
    gamma = s * sqrt((theta / s)**2 - (dx / s) * (dp / s))
    if (stp < stx)      gamma = -gamma
    p = (gamma - dx) + theta
    q = ( (gamma - dx) + gamma ) + dp
    r = p / q
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx)+dx)) / 2.0) * (stp - stx)
    if (abs(stpc - stx) < abs(stpq - stx)) then
        stpf = stpc
    else
        stpf = stpc + (stpq - stpc) / 2.0;
    end if
    this%bracketed = .true.
! Second case. A lower function value and derivatives of opposite sign. The
! minimum is bracketed. If the cubic step is closer to stx than the quadratic
! (secant) step, the cubic step is taken, else the quadratic step is taken.
else if (sgnd < 0.0) then
    info = 2
    bound = .false.
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
    s = max(abs(theta), max(abs(dx), abs(dp)))
    gamma = s * sqrt((theta / s)**2 - (dx / s) * (dp / s));
    if (stp > stx) gamma = -gamma
    p = (gamma - dp) + theta
    q = ((gamma - dp) + gamma) + dx
    r = p/q
    stpc = stp + r * (stx - stp)
    stpq = stp + (dp / (dp - dx)) * (stx - stp)
    if (abs(stpc - stp) > abs(stpq - stp)) then
        stpf = stpc
    else
        stpf = stpq
    end if
    this%bracketed = .true.
! Third case. A lower function value, derivatives of the same sign, and the
! magnitude of the derivative decreases. The cubic step is only used if the
! cubic tends to inf0y in the direction of the step or if the minimum of the
! cubic is beyond stp. Otherwise the cubic step is defined to be either
! curMinStep or curMaxStep. Thequadratic (secant) step is also computed and if
! the minimum is bracketed then thethe step closest to stx is taken, else the
! step farthest away is taken.
else if(abs(dp) <abs(dx)) then
    info = 3
    bound = .true.
    theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp
    s = max(abs(theta), max(abs(dx), abs(dp)))
    ! The case gamma = 0 only arises if the cubic does not tend to inf0y in the
    ! direction of the step.
    gamma = s*sqrt(max(0.0,(theta / s)**2 - (dx / s) * (dp / s)))
    if (stp > stx)  gamma = -gamma
    p = (gamma - dp) + theta
    q = (gamma + (dx - dp)) + gamma
    r = p/q
    if (r < 0.0 .and. gamma /= 0.0) then
        stpc = stp + r*(stx - stp)
    else if (stp > stx) then
        stpc = this%curMaxStep
    else
        stpc = this%curMinStep
    end if

    stpq = stp + (dp/(dp-dx))*(stx - stp)
    if (this%bracketed) then
        if (abs(stp-stpc) < abs(stp-stpq)) then
            stpf = stpc
        else
            stpf = stpq
        end if
    else
        if (abs(stp-stpc) > abs(stp-stpq)) then
            stpf = stpc
        else
            stpf = stpq
        end if
    end if
! Fourth case. A lower function value, derivatives of the same sign, and the
! magnitude of the derivative does not decrease. If the minimum is not
! bracketed, the step is either curMinStep or curMaxStep, else the cubic step
! is taken.
else
    info = 4
    bound = .false.
    if (this%bracketed) then
        theta = 3*(fp - fy)/(sty - stp) + dy + dp
        s = max(abs(theta), max(abs(dy),abs(dp)))
        gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
        if (stp > sty) gamma = -gamma;
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dy
        r = p/q
        stpc = stp + r*(sty - stp)
        stpf = stpc
    else if (stp > stx) then
        stpf = this%curMaxStep
    else
        stpf = this%curMinStep
    end if
end if

! Update the interval of uncertainty. This update does not depend on the new
! step or the case analysis above.
if (fp > fx) then
    sty = stp
    fy = fp
    dy = dp
else
    if (sgnd < 0.0) then
        sty = stx
        fy = fx
        dy = dx
    end if
    stx = stp
    fx = fp
    dx = dp
end if

! Compute the new step and safeguard it.
stpf = min(this%curMaxStep, stpf)
stpf = max(this%curMinStep, stpf)
stp  = stpf
if (this%bracketed .and. bound) then
    if (sty > stx) then
        stp = min(stx + (2.0/3.0) * (sty-stx), stp)
    else
        stp = max(stx + (2.0/3.0) * (sty-stx), stp)
    end if
end if

take_step_success = info > 0

end subroutine take_step

end module line_search
