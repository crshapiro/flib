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
module lbfgsb
!*******************************************************************************
! This module minimizes a function using the L-BFGS-B algorithm.
!
! This is a converted version of the original F77 source code published under
! the New BSD license by Jorge Nocedal and Jose Luis Morales. The original
! license header is shown below. Please read attached file lbfgsb_license.txt
! Requires Blas and line_search algorithm
! Some Linpack files are also included.
!
! See: http://users.iems.northwestern.edu/~nocedal/lbfgsb.html
!
!    L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”
!    or “3-clause license”)
!
!  ===========   L-BFGS-B (version 3.0.  April 25, 2011  =======================
!
!       This is a modified version of L-BFGS-B.
!
!       Major changes are described in the accompanying paper:
!
!           Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
!           L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
!           Optimization"  (2011). To appear in  ACM Transactions on
!           Mathematical Software,
!
!       The paper describes an improvement and a correction to Algorithm 778.
!       It is shown that the performance of the algorithm can be improved
!       significantly by making a relatively simple modication to the subspace
!       minimization phase. The correction concerns an error caused by the use
!       of routine dpmeps to estimate machine precision.
!
!              J. Nocedal  Department of Electrical Engineering and
!                          Computer Science.
!                          Northwestern University. Evanston, IL. USA
!
!
!             J.L Morales  Departamento de Matematicas,
!                          Instituto Tecnologico Autonomo de Mexico
!                          Mexico D.F. Mexico.
!
!                          March  2011
!
use stl
use minimize
use line_search

private

public :: lbfgsb_t
type :: lbfgsb_t
    ! Pointer to a minimize class
    class(minimize_t), pointer :: mini => NULL()
    ! maximum number of iterations
    integer :: maxiteri = 10000
    ! array of lower bounds
    real(rprec), dimension(:), allocatable :: lb
    ! array of upper bounds
    real(rprec), dimension(:), allocatable :: ub
    ! convergence level
    real(rprec) :: tol = 1E-6
    ! line search class
    type(line_search_t) :: ls
contains
    procedure, public :: minimize
    procedure, private :: minimize_priv
end type lbfgsb_t

interface lbfgsb_t
    module procedure :: constructor
end interface lbfgsb_t

contains

!*******************************************************************************
function constructor(i_mini, i_maxiteri, i_lb, i_ub, i_tol) result(this)
!*******************************************************************************
! Constructor for L-BFGS-B that takes as an argument a pointer to a minimize_t
!
type(lbfgsb_t) :: this
class(minimize_t), target :: i_mini
integer, intent(in), optional :: i_maxiteri
real(rprec), intent(in), dimension(:), optional :: i_lb, i_ub
real(rprec), intent(in), optional :: i_tol

! Assign input arguments
if ( present(i_maxiteri) ) this%maxiteri = i_maxiteri
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
this%ls = line_search_t(i_mini, .false., 1E-3_rprec, 0.9_rprec, 0.1_rprec)
end function constructor

!*******************************************************************************
subroutine minimize(this, x)
!*******************************************************************************
! Minimize the function and return the result in the array x
!
class(lbfgsb_t), intent(inout) :: this
real(rprec), dimension(:), intent(inout) :: x
integer :: n, m

! Set sizes of arrays
n = size(x)
m = 20

! Call private method
call this%minimize_priv(x, n, m)

end subroutine minimize

!*******************************************************************************
subroutine minimize_priv(this, x, n, m)
!*******************************************************************************
! Private minimization routine. This implementation uses the following
! local variables:
!
! n is the dimension of the problem.
! m is the maximum number of variable metric corrections
!   used to define the limited memory matrix.
! x is an approximation to the solution.
! l is the lower bound on x.
! u is the upper bound on x.
! nbd represents the type of bounds imposed on the
!     variables, and must be specified as follows:
!     nbd(i)=0 if x(i) is unbounded,
!            1 if x(i) has only a lower bound,
!            2 if x(i) has both lower and upper bounds, and
!            3 if x(i) has only an upper bound.
! f is the value of the function at x.
! g is the value of the gradient at x.
! factr: The iteration will stop
!     (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!     where epsmch is the machine precision, which is automatically
!     generated by the code. Typical values for factr: 1.d+12 for
!     low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
!     high accuracy.
! pgtol: The iteration will stop when
!             max{|proj g_i | i = 1, ..., n} <= pgtol
!     where pg_i is the ith component of the projected gradient.
! wa is a working array
! iwa is a working array
! task is a working string of characters of length 60 indicating
!   the current job when entering and quitting this subroutine.
! iprint is an integer variable that must be set by the user.
!   It controls the frequency and type of output generated:
!    iprint<0    no output is generated;
!    iprint=0    print only 1._rprec line at the last iteration;
!    0<iprint<99 print also f and |proj g| every iprint iterations;
!    iprint=99   print details of every iteration except n-vectors;
!    iprint=100  print also the changes of active set and final x;
!    iprint>100  print details of every iteration including x and g;
!   When iprint > 0, the file iteriate.dat will be created to
!                    summarize the iteration.
!
class(lbfgsb_t), intent(inout) :: this
integer, intent(in) :: n, m
real(rprec), dimension(n), intent(inout) :: x
real(rprec), dimension(n) :: l, u, g
! real(rprec), dimension((2*m + 5)*n + 11*m*m + 8*m) :: wa
real(rprec), dimension(n, m) :: ws, wy
real(rprec), dimension(m, m) :: sy
! integer, dimension(3*n) :: iwa
integer, dimension(n) :: nbd
real(rprec) :: f, factr, pgtol
character*60 :: task
integer :: iprint
real(rprec), dimension(8*m) :: wa
real(rprec), dimension(m, m) :: ss, wt
real(rprec), dimension(2*m, 2*m) :: wn, snd
real(rprec), dimension(n) :: z, r, d, t, xp
integer, dimension(n) :: index, iwhere, indx2
logical :: prjctd, cnstnd, boxed, updatd, wrk
character*3 :: word
integer :: i, k, nintol, iback, nskip, head, col, iter, itail, iupdat
integer :: nseg, nfgv, info, ifun, iword, nfree, nact, ileave, nenter
real(rprec) :: theta, fold, ddot, dr, rr, tol, xstep, sbgnrm, ddum, dnorm, dtd
real(rprec) :: epsmch, cpu1, cpu2, cachyt, sbtime, lnscht, time1, time2, gd
real(rprec) :: gdold, stp, stpmx, time, a1, a2
real(rprec), parameter :: big = 1.0E10_rprec
integer :: dummy

! Initialize
if (allocated(this%lb)) then
    l = this%lb
else
    l = -1000000000._rprec
end if
if (allocated(this%ub)) then
    u = this%ub
else
    u = 1000000000._rprec
end if
nbd = 2
epsmch = epsilon(1._rprec)
factr = this%tol/epsmch
pgtol = 0._rprec
task = 'START'
iprint = -1
call cpu_time(time1)

! Initialize counters and scalars when task='START'.
! for the limited memory BFGS matrices:
col    = 0
head   = 1
theta  = 1._rprec
iupdat = 0
updatd = .false.
iback  = 0
itail  = 0
iword  = 0
nact   = 0
ileave = 0
nenter = 0
fold   = 0._rprec
dnorm  = 0._rprec
cpu1   = 0._rprec
gd     = 0._rprec
stpmx  = 0._rprec
sbgnrm = 0._rprec
stp    = 0._rprec
gdold  = 0._rprec
dtd    = 0._rprec
! for operation counts:
iter   = 0
nfgv   = 0
nseg   = 0
nintol = 0
nskip  = 0
nfree  = n
ifun   = 0
! for stopping tolerance:
tol = factr*epsmch

! for measuring running time:
cachyt = 0
sbtime = 0
lnscht = 0

! 'word' records the status of subspace solutions.
word = '---'

! 'info' records the termination information.
info = 0

! Check the input arguments for errors.
call errclb(n,m,factr,l,u,nbd,task,info,k)
if (task(1:5) .eq. 'ERROR') then
    write(*,*) "lbfgsb_t%minimize: " // "routine detected an error"
    return
end if

! Initialize iwhere & project x onto the feasible set.
call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed)

! set iteration counter, evaluate function and derivative, and begin
task = 'FG_START'
iter = 0
call this%mini%eval(x, f, g)

! Compute the infinity norm of the (-) projected gradient.
call projgr(n,l,u,nbd,x,g,sbgnrm)

if (sbgnrm .le. pgtol) then
    ! terminate the algorithm.
    task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
    call cpu_time(time2)
    time = time2 - time1
endif

! Iterate
do while (iter < this%maxiteri)
    ! When (1:4)='CONV', the termination test in L-BFGS-B has been
    !   satisfied;
    if (task(1:4) == 'CONV') then
        exit
    end if

    ! When task(1:4)='ABNO', the routine has terminated abnormally
    !   without being able to satisfy the termination conditions,
    !   x contains the best approximation found,
    !   f and g contain f(x) and g(x) respectively;
    if (task(1:4) == 'ABNO') then
        write(*,*) "lbfgsb_t%minimize: routine terminated abnormally"
        exit
    end if

    ! When task(1:5)='ERROR', the routine has detected an error in the
    !   input parameters;
    if (task(1:4) == 'ERROR') then
        write(*,*) "lbfgsb_t%minimize: routine detected an error"
        exit
    end if

    ! If we're at a new x location
    if (task(1:5) .eq. 'NEW_X') then
        ! calculate and print out the quantities related to the new X.
        call cpu_time(cpu2)
        lnscht = lnscht + cpu2 - cpu1
        iter = iter + 1

        ! Compute the infinity norm of the projected (-)gradient.
        call projgr(n,l,u,nbd,x,g,sbgnrm)

        ! Test for termination.
        if (sbgnrm .le. pgtol) then
            ! terminate the algorithm.
            task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
            call cpu_time(time2)
            time = time2 - time1
            cycle
        endif

        ddum = max(abs(fold), abs(f), 1._rprec)
        if ((fold - f) .le. tol*ddum) then
            ! terminate the algorithm.
            task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
            if (iback .ge. 10) info = -5
                ! i.e., to issue a warning if iback>10 in the line search.
            call cpu_time(time2)
            time = time2 - time1
            cycle
        endif

        ! Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
        do i = 1, n
            r(i) = g(i) - r(i)
        end do
        rr = ddot(n,r,1,r,1)
        if (stp .eq. 1._rprec) then
            dr = gd - gdold
            ddum = -gdold
        else
            dr = (gd - gdold)*stp
            call dscal(n,stp,d,1)
            ddum = -gdold*stp
        endif

        if (dr .le. epsmch*ddum) then
            ! skip the L-BFGS update.
            nskip = nskip + 1
            updatd = .false.
        else

            ! Update the L-BFGS matrix.
            updatd = .true.
            iupdat = iupdat + 1

            ! Update matrices WS and WY and form the middle matrix in B.
            call matupd(n,m,ws,wy,sy,ss,d,r,itail,iupdat,col,head,theta,rr,dr, &
                stp,dtd)

            ! Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
            !   Store T in the upper triangular of the array wt;
            !   Cholesky factorize T to J*J' with
            !   J' stored in the upper triangular of wt.
            call formt(m,wt,sy,ss,col,theta,info)

            if (info .ne. 0) then
                ! nonpositive definiteness in Cholesky factorization;
                ! refresh the lbfgs memory and restart the iteration.
                info = 0
                col = 0
                head = 1
                theta = 1._rprec
                iupdat = 0
                updatd = .false.
            endif
            ! Now the inverse of the middle matrix in B is
            !   [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
            !   [ -L*D^(-1/2)   J ] [  0        J'          ]
        end if
    end if

    if (task(1:4) .eq. 'STOP') then
        if (task(7:9) .eq. 'CPU') then
            ! restore the previous iterate.
            call dcopy(n,t,1,x,1)
            call dcopy(n,r,1,g,1)
            f = fold
        endif
        call cpu_time(time2)
        time = time2 - time1
        cycle
    endif

    iword = -1
    if (.not. cnstnd .and. col .gt. 0) then
        ! skip the search for GCP.
        call dcopy(n,x,1,z,1)
        wrk = updatd
        nseg = 0
    else
        ! Compute the Generalized Cauchy Point (GCP).
        call cpu_time(cpu1)
        call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,m,wy,ws,sy,wt,theta,col,  &
            head,wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg,iprint, sbgnrm, info,&
            epsmch)
        if (info .ne. 0) then
            ! singular triangular system detected; refresh the lbfgs memory.
            info   = 0
            col    = 0
            head   = 1
            theta  = 1._rprec
            iupdat = 0
            updatd = .false.
            call cpu_time(cpu2)
            cachyt = cachyt + cpu2 - cpu1
            cycle
        endif
        call cpu_time(cpu2)
        cachyt = cachyt + cpu2 - cpu1
        nintol = nintol + nseg

        ! Count the entering and leaving variables for iter > 0;
        call freev(n,nfree,index,nenter,ileave,indx2,iwhere,wrk,updatd,cnstnd, &
            iprint,iter)
        nact = n - nfree
    endif

    ! If there are no free variables or B=theta*I, then
    ! skip the subspace minimization.
    if (.not.(nfree .eq. 0 .or. col .eq. 0)) then
        ! Subspace minimization.
        call cpu_time(cpu1)

        ! Form  the LEL^T factorization of the indefinite
        ! matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
        !               [L_a -R_z           theta*S'AA'S ]
        ! where     E = [-I  0]
        !               [ 0  I]
        if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat,          &
                      updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
        if (info .ne. 0) then
            ! nonpositive definiteness in Cholesky factorization;
            ! refresh the lbfgs memory and restart the iteration.
            info   = 0
            col    = 0
            head   = 1
            theta  = 1._rprec
            iupdat = 0
            updatd = .false.
            call cpu_time(cpu2)
            sbtime = sbtime + cpu2 - cpu1
            cycle
        endif

        ! compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) from 'cauchy').
        call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,theta,col,head,nfree,     &
            cnstnd,info)
        if (info .eq. 0) then
            ! call the direct method.
            call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy,       &
                theta, x, g, col, head, iword, wa, wn, iprint, info)
        end if
        if (info .ne. 0) then
            ! singular triangular system detected;
            ! refresh the lbfgs memory and restart the iteration.
            info   = 0
            col    = 0
            head   = 1
            theta  = 1._rprec
            iupdat = 0
            updatd = .false.
            call cpu_time(cpu2)
            sbtime = sbtime + cpu2 - cpu1
            cycle
        endif

        call cpu_time(cpu2)
        sbtime = sbtime + cpu2 - cpu1
    endif

    ! Line search and optimality tests.
    ! Generate the search direction d:=z-x.
    do i = 1, n
        d(i) = z(i) - x(i)
    end do
    call cpu_time(cpu1)

    dtd = ddot(n,d,1,d,1)
    dnorm = sqrt(dtd)

    ! Determine the maximum step length.
    stpmx = big
    if (cnstnd) then
        if (iter .eq. 0) then
            stpmx = 1._rprec
        else
            do i = 1, n
                a1 = d(i)
                if (nbd(i) .ne. 0) then
                    if (a1 .lt. 0._rprec .and. nbd(i) .le. 2) then
                        a2 = l(i) - x(i)
                        if (a2 .ge. 0._rprec) then
                            stpmx = 0._rprec
                        else if (a1*stpmx .lt. a2) then
                            stpmx = a2/a1
                        endif
                    else if (a1 .gt. 0._rprec .and. nbd(i) .ge. 2) then
                        a2 = u(i) - x(i)
                        if (a2 .le. 0._rprec) then
                            stpmx = 0._rprec
                        else if (a1*stpmx .gt. a2) then
                            stpmx = a2/a1
                        endif
                    endif
                endif
            enddo
        endif
    endif

    if (iter .eq. 0 .and. .not. boxed) then
        stp = min(1._rprec/dnorm, stpmx)
    else
        stp = 1._rprec
    endif

    call dcopy(n,x,1,t,1)
    call dcopy(n,g,1,r,1)
    fold = f
    ifun = 0
    iback = 0

    ! Do line search
    stp = 1._rprec
    this%ls%maxStep = stpmx
    call this%mini%eval(x, f, g)
    write(*,*) f
    gdold = ddot(n,g,1,d,1)
    if (sum(g*d) > 0.0) cycle
    call this%ls%search(x, f, g, d, stp, dummy)
    gd = ddot(n,g,1,d,1)
    nfgv = nfgv + dummy
    xstep = stp*dnorm
    task = 'NEW_X'

end do

! Evaluate minimization at current point
call this%mini%eval(x, f, g)

! Set output if present
! if ( present(o_x) ) o_x = this%x

! Print result
write(*,*) 'L-BFGS-B terminated after ', iter, 'iterations. Minimum f = ',f

end subroutine minimize_priv

!*******************************************************************************
subroutine dpofa(a,lda,n,info)
!*******************************************************************************
! factors a double precision symmetric positive definite matrix.
!
! dpofa is usually called by dpoco, but it can be called
! directly with a saving in time if  rcond  is not needed.
! (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!
! on entry
!
!    a       double precision(lda, n)
!            the symmetric matrix to be factored.  only the
!            diagonal and upper triangle are used.
!
!    lda     integer
!            the leading dimension of the array  a .
!
!    n       integer
!            the order of the matrix  a .
!
! on return
!
!    a       an upper triangular matrix  r  so that  a = trans(r)*r
!            where  trans(r)  is the transpose.
!            the strict lower triangle is unaltered.
!            if  info .ne. 0 , the factorization is not complete.
!
!    info    integer
!            = 0  for normal return.
!            = k  signals an error condition.  the leading minor
!                 of order  k  is not positive definite.
!
! linpack.  this version dated 08/14/78 .
! cleve moler, university of new mexico, argonne national lab.
!
integer lda,n,info
double precision a(lda,*)
double precision ddot,t
double precision s
integer j,jm1,k

do j = 1, n
info = j
s = 0.0d0
jm1 = j - 1
if (jm1 .ge. 1) then
    do k = 1, jm1
        t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
        t = t/a(k,k)
        a(k,j) = t
        s = s + t*t
    end do
end if
s = a(j,j) - s
if (s .le. 0.0d0) return
a(j,j) = sqrt(s)
end do
info = 0

end subroutine dpofa

!*******************************************************************************
subroutine dtrsl(t,ldt,n,b,job,info)
!*******************************************************************************
! solves systems of the form
!
!               t * x = b
! or
!               trans(t) * x = b
!
! where t is a triangular matrix of order n. here trans(t)
! denotes the transpose of the matrix t.
!
! on entry
!
!     t         double precision(ldt,n)
!               t contains the matrix of the system. the 0.0d0
!               elements of the matrix are not referenced, and
!               the corresponding elements of the array can be
!               used to store other information.
!
!     ldt       integer
!               ldt is the leading dimension of the array t.
!
!     n         integer
!               n is the order of the system.
!
!     b         double precision(n).
!               b contains the right hand side of the system.
!
!     job       integer
!               job specifies what kind of system is to be solved.
!               if job is
!
!                    00   solve t*x=b, t lower triangular,
!                    01   solve t*x=b, t upper triangular,
!                    10   solve trans(t)*x=b, t lower triangular,
!                    11   solve trans(t)*x=b, t upper triangular.
!
! on return
!
!     b         b contains the solution, if info .eq. 0.
!               otherwise b is unaltered.
!
!     info      integer
!               info contains 0.0d0 if the system is nonsingular.
!               otherwise info contains the index of
!               the first 0.0d0 diagonal element of t.
!
! linpack. this version dated 08/14/78 .
! g. w. stewart, university of maryland, argonne national lab.
!
integer ldt,n,job,info
double precision t(ldt,*),b(*)
double precision ddot,temp
integer c, j, jj

!heck for 0.0d0 diagonal elements.
do info = 1, n
    if (t(info,info) .eq. 0.0d0) return
end do
info = 0

! determine the task and go to it.
c = 1
if (mod(job,10) .ne. 0) c = 2
if (mod(job,100)/10 .ne. 0) c = c + 2

select case (c)

    case (1)
        ! solve t*x=b for t lower triangular
        b(1) = b(1)/t(1,1)
        if (n .lt. 2) return
        do j = 2, n
           temp = -b(j-1)
           call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
           b(j) = b(j)/t(j,j)
        end do

    case (2)
        ! solve t*x=b for t upper triangular.
        b(n) = b(n)/t(n,n)
        if (n .lt. 2) return
        do jj = 2, n
            j = n - jj + 1
            temp = -b(j+1)
            call daxpy(j,temp,t(1,j+1),1,b(1),1)
            b(j) = b(j)/t(j,j)
        end do

    case (3)
        ! solve trans(t)*x=b for t lower triangular.
        b(n) = b(n)/t(n,n)
        if (n .lt. 2) return
        do jj = 2, n
            j = n - jj + 1
            b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
            b(j) = b(j)/t(j,j)
        end do

    case (4)
        ! solve trans(t)*x=b for t upper triangular.
        b(1) = b(1)/t(1,1)
        if (n .lt. 2) return
        do j = 2, n
            b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
            b(j) = b(j)/t(j,j)
        end do
end select

end subroutine dtrsl

!*******************************************************************************
subroutine active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)
!*******************************************************************************
! This subroutine initializes iwhere and projects the initial x to
! the feasible set if necessary.
!
! iwhere is an integer array of dimension n.
!   On entry iwhere is unspecified.
!   On exit iwhere(i)=-1  if x(i) has no bounds
!                     3   if l(i)=u(i)
!                     0   otherwise.
!   In cauchy, iwhere is given finer gradations.
!
logical          prjctd, cnstnd, boxed
integer          n, iprint, nbd(n), iwhere(n)
double precision x(n), l(n), u(n)
integer          nbdd, i

nbdd = 0
prjctd = .false.
cnstnd = .false.
boxed = .true.

! Project the initial x to the easible set if necessary.
do i = 1, n
   if (nbd(i) .gt. 0) then
      if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
         if (x(i) .lt. l(i)) then
            prjctd = .true.
            x(i) = l(i)
         endif
         nbdd = nbdd + 1
      else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
         if (x(i) .gt. u(i)) then
            prjctd = .true.
            x(i) = u(i)
         endif
         nbdd = nbdd + 1
      endif
   endif
end do

! Initialize iwhere and assign values to cnstnd and boxed.
do i = 1, n
    if (nbd(i) .ne. 2) boxed = .false.
    if (nbd(i) .eq. 0) then
        ! this variable is always free
        iwhere(i) = -1
        ! otherwise set x(i)=mid(x(i), u(i), l(i)).
   else
       cnstnd = .true.
       if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. 0.0d0) then
           ! this variable is always fixed
           iwhere(i) = 3
       else
           iwhere(i) = 0
       endif
    endif
end do

if (iprint .ge. 0) then
    if (prjctd) write (6,*) "The initial X is infeasible.  Restart with its projection."
    if (.not. cnstnd) write (6,*) "This problem is unconstrained."
endif

if (iprint .gt. 0) then
    write (6,"(/,'At X0 ',i9,' variables are exactly at the bounds')") nbdd
end if

end subroutine active

!*******************************************************************************
subroutine bmv(m, sy, wt, col, v, p, info)
!*******************************************************************************
! This subroutine computes the product of the 2m x 2m middle matrix
!   in the compact L-BFGS formula of B and a 2m vector v;
!   it returns the product in p.
!
! m is an integer variable.
!   On entry m is the maximum number of variable metric corrections
!     used to define the limited memory matrix.
!   On exit m is unchanged.
!
! sy is a double precision array of dimension m x m.
!   On entry sy specifies the matrix S'Y.
!   On exit sy is unchanged.
!
! wt is a double precision array of dimension m x m.
!   On entry wt specifies the upper triangular matrix J' which is
!     the Cholesky factor of (thetaS'S+LD^(-1)L').
!   On exit wt is unchanged.
!
! col is an integer variable.
!   On entry col specifies the number of s-vectors (or y-vectors)
!     stored in the compact L-BFGS formula.
!   On exit col is unchanged.
!
! v is a double precision array of dimension 2col.
!   On entry v specifies vector v.
!   On exit v is unchanged.
!
! p is a double precision array of dimension 2col.
!   On entry p is unspecified.
!   On exit p is the product Mv.
!
! info is an integer variable.
!   On entry info is unspecified.
!   On exit info = 0       for normal return,
!                = nonzero for abnormal return when the system
!                            to be solved by dtrsl is singular.
!
integer m, col, info
double precision sy(m, m), wt(m, m), v(2*col), p(2*col)
integer :: i, k, i2
double precision :: sum

if (col .eq. 0) return

! PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
!               [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
! solve Jp2=v2+LD^(-1)v1.
p(col + 1) = v(col + 1)
do i = 2, col
    i2 = col + i
    sum = 0.0d0
    do k = 1, i - 1
        sum = sum + sy(i,k)*v(k)/sy(k,k)
    enddo
    p(i2) = v(i2) + sum
enddo
! Solve the triangular system
call dtrsl(wt,m,col,p(col+1),11,info)
if (info .ne. 0) return

! olve D^(1/2)p1=v1.
do i = 1, col
    p(i) = v(i)/sqrt(sy(i,i))
enddo

! PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
!                [  0         J'           ] [ p2 ]   [ p2 ].
! solve J^Tp2=p2.
call dtrsl(wt,m,col,p(col+1),01,info)
if (info .ne. 0) return

! compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
!           =-D^(-1/2)p1+D^(-1)L'p2.
do i = 1, col
    p(i) = -p(i)/sqrt(sy(i,i))
enddo

do i = 1, col
    sum = 0.d0
    do k = i + 1, col
        sum = sum + sy(k,i)*p(col+k)/sy(i,i)
    enddo
    p(i) = p(i) + sum
enddo

end subroutine bmv

!*******************************************************************************
subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, m, wy, ws, sy,&
    wt, theta, col, head, p, c, wbp,v, nseg, iprint, sbgnrm, info, epsmch)
!*******************************************************************************
! For given x, l, u, g (with sbgnrm > 0), and a limited memory
!   BFGS matrix B defined in terms of matrices WY, WS, WT, and
!   scalars head, col, and theta, this subroutine computes the
!   generalized Cauchy point (GCP), defined as the first local
!   minimizer of the quadratic
!
!              Q(x + s) = g's + 1/2 s'Bs
!
!   along the projected gradient direction P(x-tg,l,u).
!   The routine returns the GCP in xcp.
!
! n is an integer variable.
!   On entry n is the dimension of the problem.
!   On exit n is unchanged.
!
! x is a double precision array of dimension n.
!   On entry x is the starting point for the GCP computation.
!   On exit x is unchanged.
!
! l is a double precision array of dimension n.
!   On entry l is the lower bound of x.
!   On exit l is unchanged.
!
! u is a double precision array of dimension n.
!   On entry u is the upper bound of x.
!   On exit u is unchanged.
!
! nbd is an integer array of dimension n.
!   On entry nbd represents the type of bounds imposed on the
!     variables, and must be specified as follows:
!     nbd(i)=0 if x(i) is unbounded,
!            1 if x(i) has only a lower bound,
!            2 if x(i) has both lower and upper bounds, and
!            3 if x(i) has only an upper bound.
!   On exit nbd is unchanged.
!
! g is a double precision array of dimension n.
!   On entry g is the gradient of f(x).  g must be a nonzero vector.
!   On exit g is unchanged.
!
! iorder is an integer working array of dimension n.
!   iorder will be used to store the breakpoints in the piecewise
!   linear path and free variables encountered. On exit,
!     iorder(1),...,iorder(nleft) are indices of breakpoints
!                            which have not been encountered;
!     iorder(nleft+1),...,iorder(nbreak) are indices of
!                                 encountered breakpoints; and
!     iorder(nfree),...,iorder(n) are indices of variables which
!             have no bound constraits along the search direction.
!
! iwhere is an integer array of dimension n.
!   On entry iwhere indicates only the permanently fixed (iwhere=3)
!   or free (iwhere= -1) components of x.
!   On exit iwhere records the status of the current x variables.
!   iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
!             0   if x(i) is free and has bounds, and is moved
!             1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
!             2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
!             3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
!             -1  if x(i) is always free, i.e., it has no bounds.
!
! t is a double precision working array of dimension n.
!   t will be used to store the break points.
!
! d is a double precision array of dimension n used to store
!   the Cauchy direction P(x-tg)-x.
!
! xcp is a double precision array of dimension n used to return the
!   GCP on exit.
!
! m is an integer variable.
!   On entry m is the maximum number of variable metric corrections
!     used to define the limited memory matrix.
!   On exit m is unchanged.
!
! ws, wy, sy, and wt are double precision arrays.
!   On entry they store information that defines the
!                         limited memory BFGS matrix:
!     ws(n,m) stores S, a set of s-vectors;
!     wy(n,m) stores Y, a set of y-vectors;
!     sy(m,m) stores S'Y;
!     wt(m,m) stores the
!             Cholesky factorization of (theta*S'S+LD^(-1)L').
!   On exit these arrays are unchanged.
!
! theta is a double precision variable.
!   On entry theta is the scaling factor specifying B_0 = theta I.
!   On exit theta is unchanged.
!
! col is an integer variable.
!   On entry col is the actual number of variable metric
!     corrections stored so far.
!   On exit col is unchanged.
!
! head is an integer variable.
!   On entry head is the location of the first s-vector (or y-vector)
!     in S (or Y).
!   On exit col is unchanged.
!
! p is a double precision working array of dimension 2m.
!   p will be used to store the vector p = W^(T)d.
!
! c is a double precision working array of dimension 2m.
!   c will be used to store the vector c = W^(T)(xcp-x).
!
! wbp is a double precision working array of dimension 2m.
!   wbp will be used to store the row of W corresponding
!     to a breakpoint.
!
! v is a double precision working array of dimension 2m.
!
! nseg is an integer variable.
!   On exit nseg records the number of quadratic segments explored
!     in searching for the GCP.
!
! sg and yg are double precision arrays of dimension m.
!   On entry sg  and yg store S'g and Y'g correspondingly.
!   On exit they are unchanged.
!
! iprint is an INTEGER variable that must be set by the user.
!   It controls the frequency and type of output generated:
!    iprint<0    no output is generated;
!    iprint=0    print only 1.0d0 line at the last iteration;
!    0<iprint<99 print also f and |proj g| every iprint iterations;
!    iprint=99   print details of every iteration except n-vectors;
!    iprint=100  print also the changes of active set and final x;
!    iprint>100  print details of every iteration including x and g;
!   When iprint > 0, the file iterate.dat will be created to
!                    summarize the iteration.
!
! sbgnrm is a double precision variable.
!   On entry sbgnrm is the norm of the projected gradient at x.
!   On exit sbgnrm is unchanged.
!
! info is an integer variable.
!   On entry info is 0.
!   On exit info = 0       for normal return,
!                = nonzero for abnormal return when the the system
!                          used in routine bmv is singular.
!
integer :: n, m, head, col, nseg, iprint, info, nbd(n), iorder(n), iwhere(n)
double precision :: theta, epsmch, x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
double precision :: wy(n, col), ws(n, col), sy(m, m)
double precision :: wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)
logical :: xlower, xupper, bnded
integer :: i, j, col2, nfree, nbreak, pointr, ibp, nleft, ibkmin, iter
double precision :: f1, f2, dt, dtm, tsum, dibp, zibp, dibp2, bkmin, tu, tl
double precision :: wmc, wmp, wmw, ddot, tj, tj0, neggi, sbgnrm, f2_org
logical :: iterate = .true., finish_loop = .true.

! initialize to prevent warning at compile time
nleft = -100000
tu = -1000000d0
tl = -1000000d0

! Check the status of the variables, reset iwhere(i) if necessary;
!   compute the Cauchy direction d and the breakpoints t; initialize
!   the derivative f1 and the vector p = W'd (for theta = 1).
if (sbgnrm .le. 0.0d0) then
    if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
    call dcopy(n,x,1,xcp,1)
    return
endif
bnded = .true.
nfree = n + 1
nbreak = 0
ibkmin = 0
bkmin = 0.0d0
col2 = 2*col
f1 = 0.0d0
if (iprint .ge. 99) write(6,*) "---------------- CAUCHY entered-------------------"

! We set p to 0.0d0 and build it up as we determine d.
do i = 1, col2
    p(i) = 0.0d0
enddo

! In the following loop we determine for each variable its bound
!   status and its breakpoint, and update p accordingly.
!   Smallest breakpoint is identified.

do i = 1, n
    neggi = -g(i)
    if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
        ! if x(i) is not a constant and has bounds,
        ! compute the difference between x(i) and its bounds.
        if (nbd(i) .le. 2) tl = x(i) - l(i)
        if (nbd(i) .ge. 2) tu = u(i) - x(i)

        ! If a variable is close enough to a bound
        !   we treat it as at bound.
        xlower = nbd(i) .le. 2 .and. tl .le. 0.0d0
        xupper = nbd(i) .ge. 2 .and. tu .le. 0.0d0

        ! reset iwhere(i).
        iwhere(i) = 0
        if (xlower) then
            if (neggi .le. 0.0d0) iwhere(i) = 1
        else if (xupper) then
            if (neggi .ge. 0.0d0) iwhere(i) = 2
        else
            if (abs(neggi) .le. 0.0d0) iwhere(i) = -3
        endif
    endif
    pointr = head
    if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
        d(i) = 0.0d0
    else
        d(i) = neggi
        f1 = f1 - neggi*neggi
        ! calculate p := p - W'e_i* (g_i).
        do j = 1, col
            p(j) = p(j) +  wy(i,pointr)* neggi
            p(col + j) = p(col + j) + ws(i,pointr)*neggi
            pointr = mod(pointr,m) + 1
        enddo
        if (nbd(i) .le. 2 .and. nbd(i) .ne. 0 .and. neggi .lt. 0.0d0) then
            ! x(i) + d(i) is bounded; compute t(i).
            nbreak = nbreak + 1
            iorder(nbreak) = i
            t(nbreak) = tl/(-neggi)
            if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
            endif
        else if (nbd(i) .ge. 2 .and. neggi .gt. 0.0d0) then
            ! x(i) + d(i) is bounded; compute t(i).
            nbreak = nbreak + 1
            iorder(nbreak) = i
            t(nbreak) = tu/neggi
            if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
            endif
        else
            ! x(i) + d(i) is not bounded.
            nfree = nfree - 1
            iorder(nfree) = i
            if (abs(neggi) .gt. 0.0d0) bnded = .false.
        endif
    endif
enddo

! The indices of the nonzero components of d are now stored
!   in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
!   The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.

if (theta .ne. 1.0d0) then
    ! complete the initialization of p for theta not= 1.0d0.
    call dscal(col,theta,p(col+1),1)
endif

! Initialize GCP xcp = x.
call dcopy(n,x,1,xcp,1)
if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
    ! is a 0.0d0 vector, return with the initial xcp as GCP.
    if (iprint .gt. 100) write (6,"('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))") (xcp(i), i = 1, n)
    return
endif

! Initialize c = W'(xcp - x) = 0.
do j = 1, col2
    c(j) = 0.0d0
enddo

! Initialize derivative f2.
f2 =  -theta*f1
f2_org  =  f2
if (col .gt. 0) then
    call bmv(m,sy,wt,col,p,v,info)
    if (info .ne. 0) return
    f2 = f2 - ddot(col2,v,1,p,1)
endif
dtm = -f1/f2
tsum = 0.0d0
nseg = 1
if (iprint .ge. 99)   write (6,*) 'There are ',nbreak,'  breakpoints '

! If there are no breakpoints, locate the GCP and return.
if (nbreak .eq. 0) then
    iterate = .false.
else
    nleft = nbreak
    iter = 1
    tj = 0.0d0
endif

do while (iterate)
    ! Find the next smallest breakpoint;
    !   compute dt = t(nleft) - t(nleft + 1).
    tj0 = tj
    if (iter .eq. 1) then
        ! Since we already have the smallest breakpoint we need not do
        ! heapsort yet. Often only 1.0d0 breakpoint is used and the
        ! cost of heapsort is avoided.
        tj = bkmin
        ibp = iorder(ibkmin)
    else
        if (iter .eq. 2) then
            ! Replace the already used smallest breakpoint with the
            ! breakpoint numbered nbreak > nlast, before heapsort call.
            if (ibkmin .ne. nbreak) then
                t(ibkmin) = t(nbreak)
                iorder(ibkmin) = iorder(nbreak)
            endif
        endif
        ! Update heap structure of breakpoints
        !   (if iter=2, initialize heap).
        call hpsolb(nleft,t,iorder,iter-2)
        tj = t(nleft)
        ibp = iorder(nleft)
    endif

    dt = tj - tj0
    if (dt .ne. 0.0d0 .and. iprint .ge. 100) then
        write (6,"(/,'Piece    ',i3,' --f1, f2 at start point ',        1p,2(1x,d11.4))") nseg,f1,f2
        write (6,"('Distance to the next break point =  ',1p,d11.4)") dt
        write (6,"('Distance to the stationary point =  ',1p,d11.4)") dtm
    endif

    ! If a minimizer is within this interval, locate the GCP and return.
    if (dtm .lt. dt) exit

    ! Otherwise fix 1.0d0 variable and
    !   reset the corresponding component of d to 0.0d0.
    tsum = tsum + dt
    nleft = nleft - 1
    iter = iter + 1
    dibp = d(ibp)
    d(ibp) = 0.0d0
    if (dibp .gt. 0.0d0) then
        zibp = u(ibp) - x(ibp)
        xcp(ibp) = u(ibp)
        iwhere(ibp) = 2
    else
        zibp = l(ibp) - x(ibp)
        xcp(ibp) = l(ibp)
        iwhere(ibp) = 1
    endif
    if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
    if (nleft .eq. 0 .and. nbreak .eq. n) then
        ! all n variables are fixed, return with xcp as GCP.
        dtm = dt
        finish_loop = .false.
        exit
    endif

    ! Update the derivative information.
    nseg = nseg + 1
    dibp2 = dibp**2

    ! Update f1 and f2.
    ! temporarily set f1 and f2 for col=0.
    f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
    f2 = f2 - theta*dibp2

    if (col .gt. 0) then
        ! update c = c + dt*p.
        call daxpy(col2,dt,p,1,c,1)

        ! choose wbp, the row of W corresponding to the breakpoint encountered.
        pointr = head
        do j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
        enddo

        ! compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
        call bmv(m,sy,wt,col,wbp,v,info)
        if (info .ne. 0) return
        wmc = ddot(col2,c,1,v,1)
        wmp = ddot(col2,p,1,v,1)
        wmw = ddot(col2,wbp,1,v,1)

        ! update p = p - dibp*wbp.
        call daxpy(col2,-dibp,wbp,1,p,1)

        ! complete updating f1 and f2 while col > 0.
        f1 = f1 + dibp*wmc
        f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
    endif

    f2 = max(epsmch*f2_org,f2)
    if (nleft .gt. 0) then
        dtm = -f1/f2
        cycle
        ! to repeat the loop for unsearched intervals.
    else if(bnded) then
        f1 = 0.0d0
        f2 = 0.0d0
        dtm = 0.0d0
    else
        dtm = -f1/f2
    endif

    if (nleft .le. 0) iterate = .false.
end do

if (finish_loop) then
    if (iprint .ge. 99) then
        write (6,*)
        write (6,*) 'GCP found in this segment'
        write (6,"('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))") nseg,f1,f2
        write (6,"('Distance to the stationary point =  ',1p,d11.4)") dtm
    endif
    if (dtm .le. 0.0d0) dtm = 0.0d0
    tsum = tsum + dtm

    ! Move free variables (i.e., the ones w/o breakpoints) and
    !   the variables whose breakpoints haven't been reached.
    call daxpy(n,tsum,d,1,xcp,1)
endif

! Update c = c + dtm*p = W'(x^c - x)
!   which will be used in computing r = Z'(B(x^c - x) + g).
if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
if (iprint .gt. 100) write (6,"('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))") (xcp(i),i = 1,n)
if (iprint .ge. 99) write (6,*) "---------------- exit CAUCHY----------------------"

end subroutine cauchy

!*******************************************************************************
subroutine formt(m, wt, sy, ss, col, theta, info)
!*******************************************************************************
!   This subroutine forms the upper half of the pos. def. and symm.
!     T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
!     of the array wt, and performs the Cholesky factorization of T
!     to produce J*J', with J' stored in the upper triangle of wt.
!
integer :: m, col, info
double precision :: theta, wt(m, m), sy(m, m), ss(m, m)
integer :: i, j, k, k1
double precision ddum

! Form the upper half of  T = theta*SS + L*D^(-1)*L',
!   store T in the upper triangle of the array wt.
do j = 1, col
    wt(1,j) = theta*ss(1,j)
enddo
do i = 2, col
    do j = i, col
        k1 = min(i,j) - 1
        ddum  = 0.0d0
        do k = 1, k1
            ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
        enddo
        wt(i,j) = ddum + theta*ss(i,j)
    enddo
enddo

! Cholesky factorize T to J*J' with
!   J' stored in the upper triangle of wt.
call dpofa(wt,m,col,info)
if (info .ne. 0) then
    info = -3
endif

end subroutine formt

!*******************************************************************************
subroutine freev(n, nfree, index, nenter, ileave, indx2,iwhere, wrk, updatd,   &
    cnstnd, iprint, iter)
!*******************************************************************************
! This subroutine counts the entering and leaving variables when
! iter > 0, and finds the index set of free and active variables
! at the GCP.
!
! cnstnd is a logical variable indicating whether bounds are present
!
! index is an integer array of dimension n
! for i=1,...,nfree, index(i) are the indices of free variables
! for i=nfree+1,...,n, index(i) are the indices of bound variables
! On entry after the first iteration, index gives
!   the free variables at the previous iteration.
! On exit it gives the free variables based on the determination
!   in cauchy using the array iwhere.
!
! indx2 is an integer array of dimension n
! On entry indx2 is unspecified.
! On exit with iter>0, indx2 indicates which variables
!    have changed status since the previous iteration.
! For i= 1,...,nenter, indx2(i) have changed from bound to free.
! For i= ileave+1,...,n, indx2(i) have changed from free to bound.
!
integer :: n, nfree, nenter, ileave, iprint, iter, index(n), indx2(n), iwhere(n)
logical :: wrk, updatd, cnstnd
integer :: iact,i,k

nenter = 0
ileave = n + 1
if (iter .gt. 0 .and. cnstnd) then
    ! count the entering and leaving variables.
    do i = 1, nfree
        k = index(i)
        if (iwhere(k) .gt. 0) then
            ileave = ileave - 1
            indx2(ileave) = k
            if (iprint .ge. 100) write (6,*) 'Variable ',k,' leaves the set of free variables'
        endif
    enddo
    do i = 1 + nfree, n
        k = index(i)
        if (iwhere(k) .le. 0) then
            nenter = nenter + 1
            indx2(nenter) = k
            if (iprint .ge. 100) write (6,*) 'Variable ',k,' enters the set of free variables'
        endif
    enddo

    if (iprint .ge. 99) write (6,*) n+1-ileave,' variables leave; ',nenter,' variables enter'
endif

wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd

! Find the index set of free and active variables at the GCP.
nfree = 0
iact = n + 1
do i = 1, n
    if (iwhere(i) .le. 0) then
        nfree = nfree + 1
        index(nfree) = i
    else
        iact = iact - 1
        index(iact) = i
    endif
enddo

if (iprint .ge. 99) write (6,*) nfree,' variables are free at GCP ',iter + 1

end subroutine freev

!*******************************************************************************
subroutine hpsolb(n, t, iorder, iheap)
!*******************************************************************************
! This subroutine sorts out the least element of t, and puts the
!   remaining elements of t in a heap.
!
! n is an integer variable.
!   On entry n is the dimension of the arrays t and iorder.
!   On exit n is unchanged.
!
! t is a double precision array of dimension n.
!   On entry t stores the elements to be sorted,
!   On exit t(n) stores the least elements of t, and t(1) to t(n-1)
!     stores the remaining elements in the form of a heap.
!
! iorder is an integer array of dimension n.
!   On entry iorder(i) is the index of t(i).
!   On exit iorder(i) is still the index of t(i), but iorder may be
!     permuted in accordance with t.
!
! iheap is an integer variable specifying the task.
!   On entry iheap should be set as follows:
!     iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
!     iheap .ne. 0 if otherwise.
!   On exit iheap is unchanged.
!
! References:
!   Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
!
integer :: iheap, n, iorder(n)
double precision :: t(n)
integer :: i,j,k,indxin,indxou
double precision :: ddum,out

if (iheap .eq. 0) then
    ! Rearrange the elements t(1) to t(n) to form a heap.
    do k = 2, n
        ddum  = t(k)
        indxin = iorder(k)

        ! Add ddum to the heap.
        i = k
        do
            if (i.gt.1) then
                j = i/2
                if (ddum .lt. t(j)) then
                    t(i) = t(j)
                    iorder(i) = iorder(j)
                    i = j
                    cycle
                endif
            endif
            t(i) = ddum
            iorder(i) = indxin
            exit
        enddo
    enddo
endif

! Assign to 'out' the value of t(1), the least member of the heap,
!   and rearrange the remaining members to form a heap as
!   elements 1 to n-1 of t.
if (n .gt. 1) then
    i = 1
    out = t(1)
    indxou = iorder(1)
    ddum  = t(n)
    indxin  = iorder(n)

    ! Restore the heap
    do
        j = i+i
        if (j .le. n-1) then
            if (t(j+1) .lt. t(j)) j = j+1
            if (t(j) .lt. ddum ) then
                t(i) = t(j)
                iorder(i) = iorder(j)
                i = j
                cycle
            endif
        endif
        t(i) = ddum
        iorder(i) = indxin
        exit
    enddo

    ! Put the least member in t(n).
    t(n) = out
    iorder(n) = indxou
endif

end subroutine hpsolb

!*******************************************************************************
subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
!*******************************************************************************
! This subroutine checks the validity of the input data.
!
character*60 :: task
integer :: n, m, info, k, nbd(n)
double precision :: factr, l(n), u(n)
integer :: i

! Check the input arguments for errors.
if (n .le. 0) task = 'ERROR: N .LE. 0'
if (m .le. 0) task = 'ERROR: M .LE. 0'
if (factr .lt. 0.0d0) task = 'ERROR: FACTR .LT. 0'

! Check the validity of the arrays nbd(i), u(i), and l(i).
do i = 1, n
    if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
        ! return
        task = 'ERROR: INVALID NBD'
        info = -6
        k = i
    endif
    if (nbd(i) .eq. 2) then
        if (l(i) .gt. u(i)) then
            ! return
            task = 'ERROR: NO FEASIBLE SOLUTION'
            info = -7
            k = i
        endif
    endif
enddo

end subroutine errclb

!*******************************************************************************
subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index,                 &
    theta, col, head, nfree, cnstnd, info)
!*******************************************************************************
!   This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
!     wa(2m+1)=W'(xcp-x) from subroutine cauchy.
!
logical :: cnstnd
integer :: n, m, col, head, nfree, info, index(n)
double precision :: theta, x(n), g(n), z(n), r(n), wa(4*m)
double precision :: ws(n, m), wy(n, m), sy(m, m), wt(m, m)
integer :: i,j,k,pointr
double precision :: a1,a2

if (.not. cnstnd .and. col .gt. 0) then
    do i = 1, n
        r(i) = -g(i)
    enddo
else
    do i = 1, nfree
        k = index(i)
        r(i) = -theta*(z(k) - x(k)) - g(k)
    enddo
    call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
    if (info .ne. 0) then
        info = -8
        return
    endif
    pointr = head
    do j = 1, col
        a1 = wa(j)
        a2 = theta*wa(col + j)
        do i = 1, nfree
            k = index(i)
            r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
        enddo
        pointr = mod(pointr,m) + 1
    enddo
endif

end subroutine cmprlb

!*******************************************************************************
subroutine projgr(n, l, u, nbd, x, g, sbgnrm)
!*******************************************************************************
! This subroutine computes the infinity norm of the projected
! gradient.
!
integer :: n, nbd(n)
double precision :: sbgnrm, x(n), l(n), u(n), g(n)
integer :: i
double precision :: gi

sbgnrm = 0.0d0
do i = 1, n
    gi = g(i)
    if (nbd(i) .ne. 0) then
        if (gi .lt. 0.0d0) then
            if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
        else
            if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
        endif
    endif
    sbgnrm = max(sbgnrm,abs(gi))
enddo

end subroutine projgr

!*******************************************************************************
subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head, theta, &
    rr, dr, stp, dtd)
!*******************************************************************************
!   This subroutine updates matrices WS and WY, and forms the
!     middle matrix in B.
!
integer          n, m, itail, iupdat, col, head
double precision theta, rr, dr, stp, dtd, d(n), r(n)
double precision ws(n, m), wy(n, m), sy(m, m), ss(m, m)
integer          j, pointr
double precision ddot

! Set pointers for matrices WS and WY.
if (iupdat .le. m) then
    col = iupdat
    itail = mod(head+iupdat-2,m) + 1
else
    itail = mod(itail,m) + 1
    head = mod(head,m) + 1
endif

! Update matrices WS and WY.
call dcopy(n,d,1,ws(1,itail),1)
call dcopy(n,r,1,wy(1,itail),1)

! Set theta=yy/ys.
theta = rr/dr

! Form the middle matrix in B.
!   update the upper triangle of SS,
!   and the lower triangle of SY:
if (iupdat .gt. m) then
    ! move old information
    do j = 1, col - 1
        call dcopy(j,ss(2,j+1),1,ss(1,j),1)
        call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
    enddo
endif

! add new information: the last row of SY
! and the last column of SS:
pointr = head
do j = 1, col - 1
    sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
    ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
    pointr = mod(pointr,m) + 1
enddo
if (stp .eq. 1.0d0) then
    ss(col,col) = dtd
else
    ss(col,col) = stp*stp*dtd
endif
sy(col,col) = dr

end subroutine matupd

!*******************************************************************************
subroutine subsm (n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, theta, xx, gg, &
    col, head, iword, wv, wn, iprint, info )
!*******************************************************************************
! Given xcp, l, u, r, an index set that specifies
!   the active set at xcp, and an l-BFGS matrix B
!   (in terms of WY, WS, SY, WT, head, col, and theta),
!   this subroutine computes an approximate solution
!   of the subspace problem
!
!   (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
!
!         subject to l<=x<=u
!                   x_i=xcp_i for all i in A(xcp)
!
!   along the subspace unconstrained Newton direction
!
!      d = -(Z'BZ)^(-1) r.
!
!   The formula for the Newton direction, given the L-BFGS matrix
!   and the Sherman-Morrison formula, is
!
!      d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
!
!   where
!             K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                 [L_a -R_z           theta*S'AA'S ]
!
! Note that this procedure for computing d differs
! from that described in [1]. 1.0d0 can show that the matrix K is
! equal to the matrix M^[-1]N in that paper.
!
! n is an integer variable.
!   On entry n is the dimension of the problem.
!   On exit n is unchanged.
!
! m is an integer variable.
!   On entry m is the maximum number of variable metric corrections
!     used to define the limited memory matrix.
!   On exit m is unchanged.
!
! nsub is an integer variable.
!   On entry nsub is the number of free variables.
!   On exit nsub is unchanged.
!
! ind is an integer array of dimension nsub.
!   On entry ind specifies the coordinate indices of free variables.
!   On exit ind is unchanged.
!
! l is a double precision array of dimension n.
!   On entry l is the lower bound of x.
!   On exit l is unchanged.
!
! u is a double precision array of dimension n.
!   On entry u is the upper bound of x.
!   On exit u is unchanged.
!
! nbd is a integer array of dimension n.
!   On entry nbd represents the type of bounds imposed on the
!     variables, and must be specified as follows:
!     nbd(i)=0 if x(i) is unbounded,
!            1 if x(i) has only a lower bound,
!            2 if x(i) has both lower and upper bounds, and
!            3 if x(i) has only an upper bound.
!   On exit nbd is unchanged.
!
! x is a double precision array of dimension n.
!   On entry x specifies the Cauchy point xcp.
!   On exit x(i) is the minimizer of Q over the subspace of
!                                                    free variables.
!
! d is a double precision array of dimension n.
!   On entry d is the reduced gradient of Q at xcp.
!   On exit d is the Newton direction of Q.
!
! c    xp is a double precision array of dimension n.
!   used to safeguard the projected Newton direction
!
! c    xx is a double precision array of dimension n
!   On entry it holds the current iterate
!   On output it is unchanged
!
! c    gg is a double precision array of dimension n
!   On entry it holds the gradient at the current iterate
!   On output it is unchanged
!
! ws and wy are double precision arrays;
! theta is a double precision variable;
! col is an integer variable;
! head is an integer variable.
!   On entry they store the information defining the
!                                      limited memory BFGS matrix:
!     ws(n,m) stores S, a set of s-vectors;
!     wy(n,m) stores Y, a set of y-vectors;
!     theta is the scaling factor specifying B_0 = theta I;
!     col is the number of variable metric corrections stored;
!     head is the location of the 1st s- (or y-) vector in S (or Y).
!   On exit they are unchanged.
!
! iword is an integer variable.
!   On entry iword is unspecified.
!   On exit iword specifies the status of the subspace solution.
!     iword = 0 if the solution is in the box,
!             1 if some bound is encountered.
!
! wv is a double precision working array of dimension 2m.
!
! wn is a double precision array of dimension 2m x 2m.
!   On entry the upper triangle of wn stores the LEL^T factorization
!     of the indefinite matrix
!
!          K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!              [L_a -R_z           theta*S'AA'S ]
!                                                where E = [-I  0]
!                                                          [ 0  I]
!   On exit wn is unchanged.
!
! iprint is an INTEGER variable that must be set by the user.
!   It controls the frequency and type of output generated:
!    iprint<0    no output is generated;
!    iprint=0    print only 1.0d0 line at the last iteration;
!    0<iprint<99 print also f and |proj g| every iprint iterations;
!    iprint=99   print details of every iteration except n-vectors;
!    iprint=100  print also the changes of active set and final x;
!    iprint>100  print details of every iteration including x and g;
!   When iprint > 0, the file iterate.dat will be created to
!                    summarize the iteration.
!
! info is an integer variable.
!   On entry info is unspecified.
!   On exit info = 0       for normal return,
!                = nonzero for abnormal return
!                              when the matrix K is ill-conditioned.
!
integer          n, m, nsub, col, head, iword, iprint, info
integer          ind(nsub), nbd(n)
double precision theta
double precision l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n)
double precision ws(n, m), wy(n, m)
double precision wv(2*m), wn(2*m, 2*m)
integer          pointr, m2, col2, ibd, jy, js, i, j, k
double precision alpha, xk, dk, temp1, temp2
double precision dd_p

if (nsub .le. 0) return
if (iprint .ge. 99) write (6,*) "----------------enter SUBSM --------------------"

! Compute wv = W'Zd.
pointr = head
do i = 1, col
    temp1 = 0.0d0
    temp2 = 0.0d0
    do  j = 1, nsub
        k = ind(j)
        temp1 = temp1 + wy(k,pointr)*d(j)
        temp2 = temp2 + ws(k,pointr)*d(j)
    enddo
    wv(i) = temp1
    wv(col + i) = theta*temp2
    pointr = mod(pointr,m) + 1
enddo

! Compute wv:=K^(-1)wv.
m2 = 2*m
col2 = 2*col
call dtrsl(wn,m2,col2,wv,11,info)
if (info .ne. 0) return
do i = 1, col
    wv(i) = -wv(i)
enddo
call dtrsl(wn,m2,col2,wv,01,info)
if (info .ne. 0) return

! Compute d = (1/theta)d + (1/theta**2)Z'W wv.
pointr = head
do jy = 1, col
    js = col + jy
    do i = 1, nsub
        k = ind(i)
        d(i) = d(i) + wy(k,pointr)*wv(jy)/theta + ws(k,pointr)*wv(js)
    enddo
    pointr = mod(pointr,m) + 1
enddo

call dscal( nsub, 1.0d0/theta, d, 1 )

! Let us try the projection, d is the Newton direction
iword = 0

call dcopy ( n, x, 1, xp, 1 )
do i = 1, nsub
    k  = ind(i)
    dk = d(i)
    xk = x(k)
    if ( nbd(k) .ne. 0 ) then
        if ( nbd(k).eq.1 ) then          ! lower bounds only
            x(k) = max( l(k), xk + dk )
            if ( x(k).eq.l(k) ) iword = 1
        else
            if ( nbd(k).eq.2 ) then       ! upper and lower bounds
                xk   = max( l(k), xk + dk )
                x(k) = min( u(k), xk )
                if ( x(k).eq.l(k) .or. x(k).eq.u(k) ) iword = 1
            else
                if ( nbd(k).eq.3 ) then    ! upper bounds only
                    x(k) = min( u(k), xk + dk )
                    if ( x(k).eq.u(k) ) iword = 1
                end if
            end if
        end if
    else                                ! free variables
        x(k) = xk + dk
    end if
enddo

if ( iword.eq.0 ) then
    if (iprint .ge. 99) write (6,*) "----------------exit SUBSM --------------------"
    return
endif

! check sign of the directional derivative
dd_p = 0.0d0
do i=1, n
    dd_p  = dd_p + (x(i) - xx(i))*gg(i)
enddo
if ( dd_p .gt.0.0d0 ) then
    call dcopy( n, xp, 1, x, 1 )
!     write(6,*) ' Positive dir derivative in projection '
!     write(6,*) ' Using the backtracking step '
else
    if (iprint .ge. 99) write (6,*) "----------------exit SUBSM --------------------"
    return
endif

alpha = 1.0d0
temp1 = alpha
ibd   = 0
do i = 1, nsub
    k = ind(i)
    dk = d(i)
    if (nbd(k) .ne. 0) then
        if (dk .lt. 0.0d0 .and. nbd(k) .le. 2) then
            temp2 = l(k) - x(k)
            if (temp2 .ge. 0.0d0) then
                temp1 = 0.0d0
            else if (dk*alpha .lt. temp2) then
                temp1 = temp2/dk
            endif
        else if (dk .gt. 0.0d0 .and. nbd(k) .ge. 2) then
            temp2 = u(k) - x(k)
            if (temp2 .le. 0.0d0) then
                temp1 = 0.0d0
            else if (dk*alpha .gt. temp2) then
                temp1 = temp2/dk
            endif
        endif
        if (temp1 .lt. alpha) then
            alpha = temp1
            ibd = i
        endif
    endif
enddo

if (alpha .lt. 1.0d0) then
    dk = d(ibd)
    k = ind(ibd)
    if (dk .gt. 0.0d0) then
        x(k) = u(k)
        d(ibd) = 0.0d0
    else if (dk .lt. 0.0d0) then
        x(k) = l(k)
        d(ibd) = 0.0d0
    endif
endif
do i = 1, nsub
    k    = ind(i)
    x(k) = x(k) + alpha*d(i)
enddo

if (iprint .ge. 99) write (6,*) "----------------exit SUBSM --------------------"

end subroutine subsm

!*******************************************************************************
subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, updatd, wn, wn1, &
    m, ws, wy, sy, theta, col, head, info)
!*******************************************************************************
! This subroutine forms  the LEL^T factorization of the indefinite
!
!   matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                 [L_a -R_z           theta*S'AA'S ]
!                                                where E = [-I  0]
!                                                          [ 0  I]
! The matrix K can be shown to be equal to the matrix M^[-1]N
!   occurring in section 5.1 of [1], as well as to the matrix
!   Mbar^[-1] Nbar in section 5.3.
!
! n is an integer variable.
!   On entry n is the dimension of the problem.
!   On exit n is unchanged.
!
! nsub is an integer variable
!   On entry nsub is the number of subspace variables in free set.
!   On exit nsub is not changed.
!
! ind is an integer array of dimension nsub.
!   On entry ind specifies the indices of subspace variables.
!   On exit ind is unchanged.
!
! nenter is an integer variable.
!   On entry nenter is the number of variables entering the
!     free set.
!   On exit nenter is unchanged.
!
! ileave is an integer variable.
!   On entry indx2(ileave),...,indx2(n) are the variables leaving
!     the free set.
!   On exit ileave is unchanged.
!
! indx2 is an integer array of dimension n.
!   On entry indx2(1),...,indx2(nenter) are the variables entering
!     the free set, while indx2(ileave),...,indx2(n) are the
!     variables leaving the free set.
!   On exit indx2 is unchanged.
!
! iupdat is an integer variable.
!   On entry iupdat is the total number of BFGS updates made so far.
!   On exit iupdat is unchanged.
!
! updatd is a logical variable.
!   On entry 'updatd' is true if the L-BFGS matrix is updatd.
!   On exit 'updatd' is unchanged.
!
! wn is a double precision array of dimension 2m x 2m.
!   On entry wn is unspecified.
!   On exit the upper triangle of wn stores the LEL^T factorization
!     of the 2*col x 2*col indefinite matrix
!                 [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                 [L_a -R_z           theta*S'AA'S ]
!
! wn1 is a double precision array of dimension 2m x 2m.
!   On entry wn1 stores the lower triangular part of
!                 [Y' ZZ'Y   L_a'+R_z']
!                 [L_a+R_z   S'AA'S   ]
!     in the previous iteration.
!   On exit wn1 stores the corresponding updated matrices.
!   The purpose of wn1 is just to store these inner products
!   so they can be easily updated and inserted into wn.
!
! m is an integer variable.
!   On entry m is the maximum number of variable metric corrections
!     used to define the limited memory matrix.
!   On exit m is unchanged.
!
! ws, wy, sy, and wtyy are double precision arrays;
! theta is a double precision variable;
! col is an integer variable;
! head is an integer variable.
!   On entry they store the information defining the
!                                      limited memory BFGS matrix:
!     ws(n,m) stores S, a set of s-vectors;
!     wy(n,m) stores Y, a set of y-vectors;
!     sy(m,m) stores S'Y;
!     wtyy(m,m) stores the Cholesky factorization
!                               of (theta*S'S+LD^(-1)L')
!     theta is the scaling factor specifying B_0 = theta I;
!     col is the number of variable metric corrections stored;
!     head is the location of the 1st s- (or y-) vector in S (or Y).
!   On exit they are unchanged.
!
! info is an integer variable.
!   On entry info is unspecified.
!   On exit info =  0 for normal return;
!                = -1 when the 1st Cholesky factorization failed;
!                = -2 when the 2st Cholesky factorization failed.
!
integer :: n, nsub, m, col, head, nenter, ileave, iupdat
integer :: info, ind(n), indx2(n)
double precision :: theta, wn(2*m, 2*m), wn1(2*m, 2*m)
double precision :: ws(n, m), wy(n, m), sy(m, m)
logical :: updatd
integer :: m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k
integer :: col2,pbegin,pend,dbegin,dend,upcl
double precision :: ddot,temp1,temp2,temp3,temp4

! Form the lower triangular part of
!           WN1 = [Y' ZZ'Y   L_a'+R_z']
!                 [L_a+R_z   S'AA'S   ]
!    where L_a is the strictly lower triangular part of S'AA'Y
!          R_z is the upper triangular part of S'ZZ'Y.

if (updatd) then
    if (iupdat .gt. m) then
        ! shift old part of WN1.
        do jy = 1, m - 1
            js = m + jy
            call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
            call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
            call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
        enddo
    endif

    ! put new rows in blocks (1,1), (2,1) and (2,2).
    pbegin = 1
    pend = nsub
    dbegin = nsub + 1
    dend = n
    iy = col
    is = m + col
    ipntr = head + col - 1
    if (ipntr .gt. m) ipntr = ipntr - m
    jpntr = head
    do jy = 1, col
        js = m + jy
        temp1 = 0.0d0
        temp2 = 0.0d0
        temp3 = 0.0d0
        ! compute element jy of row 'col' of Y'ZZ'Y
        do k = pbegin, pend
            k1 = ind(k)
            temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
        enddo
        ! compute elements jy of row 'col' of L_a and S'AA'S
        do k = dbegin, dend
            k1 = ind(k)
            temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
            temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
        end do
        wn1(iy,jy) = temp1
        wn1(is,js) = temp2
        wn1(is,jy) = temp3
        jpntr = mod(jpntr,m) + 1
    end do

    ! put new column in block (2,1).
    jy = col
    jpntr = head + col - 1
    if (jpntr .gt. m) jpntr = jpntr - m
    ipntr = head
    do i = 1, col
        is = m + i
        temp3 = 0.0d0
        ! compute element i of column 'col' of R_z
        do k = pbegin, pend
            k1 = ind(k)
            temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
        enddo
        ipntr = mod(ipntr,m) + 1
        wn1(is,jy) = temp3
    enddo
    upcl = col - 1
else
    upcl = col
endif

!   modify the old parts in blocks (1,1) and (2,2) due to changes
!   in the set of free variables.
ipntr = head
do iy = 1, upcl
    is = m + iy
    jpntr = head
    do jy = 1, iy
        js = m + jy
        temp1 = 0.0d0
        temp2 = 0.0d0
        temp3 = 0.0d0
        temp4 = 0.0d0
        do k = 1, nenter
            k1 = indx2(k)
            temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
            temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
        enddo
        do k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
        end do
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
            wn1(is,js) = wn1(is,js) - temp2 + temp4
            jpntr = mod(jpntr,m) + 1
    enddo
    ipntr = mod(ipntr,m) + 1
enddo

! modify the old parts in block (2,1).
ipntr = head
do is = m + 1, m + upcl
    jpntr = head
    do jy = 1, upcl
        temp1 = 0.0d0
        temp3 = 0.0d0
        do k = 1, nenter
            k1 = indx2(k)
            temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
        end do
        do k = ileave, n
            k1 = indx2(k)
            temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
        end do
        if (is .le. jy + m) then
            wn1(is,jy) = wn1(is,jy) + temp1 - temp3
        else
            wn1(is,jy) = wn1(is,jy) - temp1 + temp3
        endif
        jpntr = mod(jpntr,m) + 1
    enddo
    ipntr = mod(ipntr,m) + 1
end do

! Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
!                                 [-L_a +R_z        S'AA'S*theta]
m2 = 2*m
do iy = 1, col
    is = col + iy
    is1 = m + iy
    do jy = 1, iy
        js = col + jy
        js1 = m + jy
        wn(jy,iy) = wn1(iy,jy)/theta
        wn(js,is) = wn1(is1,js1)*theta
    end do
    do jy = 1, iy - 1
        wn(jy,is) = -wn1(is1,jy)
    end do
    do jy = iy, col
        wn(jy,is) = wn1(is1,jy)
    end do
    wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
enddo

! Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
!                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
!
!    first Cholesky factor (1,1) block of wn to get LL'
!                      with L' stored in the upper triangle of wn.
call dpofa(wn,m2,col,info)
if (info .ne. 0) then
    info = -1
    return
endif

! then form L^-1(-L_a'+R_z') in the (1,2) block.
col2 = 2*col
do js = col+1 ,col2
    call dtrsl(wn,m2,col,wn(1,js),11,info)
end do

! Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
!    upper triangle of (2,2) block of wn.
do is = col+1, col2
    do js = is, col2
        wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
    enddo
enddo

! Cholesky factorization of (2,2) block of wn.
call dpofa(wn(col+1,col+1),m2,col,info)
if (info .ne. 0) then
    info = -2
    return
endif

end subroutine formk

end module lbfgsb
