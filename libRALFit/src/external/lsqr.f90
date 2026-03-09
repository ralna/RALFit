! July 2025
! Reverse communication variant of LSQR.
! Optionally allows split or left preconditioning of the normal equations.
! The default stopping is to use the stopping criteria of Papez and Tichy.
! The user can test for convergence 
! (or terminate after a chosen number of iterations).

! Also, allows for (partial or selective) one-sided reorthogonalization.

! Here, we base notation on the paper
! Preconditioning of LSQR and CGLS: Variants, Properties and Relations.
! Havelková E, Hnětynková I. 
! In: European Conference on Numerical Mathematics and Advanced Applications 
! (2023), pp. 415-424). Springer Nature Switzerland.
! We implement RP-LSQR, Modified RP-LSQR and LP-LSQR.

! Stopping criteria of Papez and Tichy is described in
! Estimating error norms in CG-like algorithms for least-squares and 
! least-norm problems. Jan Papez · Petr Tichy. Numerical Algorithms
! (2024) 97:1--28.

! Interface designed to be as for CGLS code, but with inclusion of
! an optional damping parameter (regularization). The action parameter
! has the same meaning so easy for user to call LSQR or CGLS
! using same driver.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! LSQR solves
!   min_x || b - Ax ||
!
!   min_x ||(   A    )*x - ( b ) ||  (regularized LS)
!         ||( damp*I )     ( 0 ) ||

! The parameter damp is intended to help regularize
! ill-conditioned systems, by preventing the true solution from
! being very large.

! If M = PP^T then right preconditioned LSQR solves
!    min_y ||b - AP^{-1} xhat|| 
!

! Right preconditioning corresponds to symmetric (split) preconditioning
! of the normal equations
!
!                P^{-1} A^T A P^{-T} xhat = P^{-1} A^T b,   P^T x_{LS} = xhat

! One choice is P = L, an incomplete factor of the normal matrix A^T A.
! This can be computed using HSL_MI35.
!
! The code also offers factorization-free preconditioned LSQR. This solves
! the LS problem by solving the left preconditioned normal equations
!
!                M^{-1} A^T A x_{LS} = M^{-1} A^T b 

! where  M \approx A^T A. Again, M = LL^T can be used but the factored form of
! M is NOT necessary.


! Note that x is not an input parameter.
    ! If some initial estimate x0 is known and if damp = 0,
    ! one could proceed as follows:
    !
    ! 1. Compute a residual vector     r0 = b - A*x0.
    ! 2. Use LSQR to solve the system  A*dx = r0.
    ! 3. Add the correction dx to obtain a final solution x = x0 + dx.
    !
    ! This requires that x0 be available before the first call 
    ! to LSQR and after the final call. 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "preprocessor.FPP"

module lsqr_reverse_double

   implicit none

#ifdef SINGLE_PRECISION
 integer, parameter :: wp = selected_real_kind(6) ! c_float
#else
 integer, parameter :: wp = selected_real_kind(15) ! c_double
#endif

   private
   public :: lsqr_keep, lsqr_options, lsqr_inform
   public :: lsqr, lsqr_free

   integer(4),  parameter :: ip = kind( 0 )
   integer, parameter :: long = selected_int_kind(18)
   real(wp),    parameter :: zero = 0.0_wp, one = 1.0_wp

   ! error flags(explained below)
   integer(ip), parameter :: lsqr_stop_m_oor           = -1
   integer(ip), parameter :: lsqr_stop_allocation      = -2
   integer(ip), parameter :: lsqr_stop_deallocation    = -3
   integer(ip), parameter :: lsqr_stop_itnlim          = -4

   ! warning flags(explained below)
   integer(ip), parameter :: lsqr_stop_anorm           = 1
   integer(ip), parameter :: lsqr_stop_x0              = 2
   integer(ip), parameter :: lsqr_stop_stagnate        = 4
   ! positive values are added if more than one warning issued

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Make interfaces generic.
   interface lsqr
      module procedure lsqr_double
   end interface

   interface lsqr_free
      module procedure lsqr_free_double  ! deallocates components of keep
   end interface

! ------------------------------------------------------------------
    
    ! The matrices A, P and M are treated as linear operators that are 
    ! accessed by reverse communication, that is, requests for matrix-vector 
    ! products with A, P or M are passed back to the user.
    !
    ! LSQR uses an iterative method to approximate the solution.
    ! The number of iterations required to reach a certain accuracy
    ! depends strongly on the scaling of the problem and on the preconditioner.

    ! Unless better information is known, we recommend that
    ! the columns of A are prescaled so that they all have
    ! the same Euclidean norm (e.g., 1.0).
    !

! ------------------------------------------------------------------
  ! control data type
  ! These may be set before first call and must not be altered between calls.
  type :: lsqr_options

    integer(ip) :: itnlim = -1 ! max number of iterations.
    ! This must be set because Papez and Tichy stopping requires arrays of
    ! size the number of iterations.
    ! If it is negative then we will use itnlim = n.

    integer(ip) :: test_frequency = 1
    ! Return to user for convergence testing every test_frequency iterations.
    ! Values less than 1 are treated as 1
    
    integer(ip) :: localSize = 0 
    ! No. of vectors for local reorthogonalization. Not used
    ! if options%orthog_choice = 0
    ! 0       No reorthogonalization is performed.
    ! >0      This many n-vectors "v" are saved for reorthogonalizing.
    !         If localSize = min(m,n) then full one-sided reorth'n is done.
    !         min(options%localSize,m,n) vectors will be allocated or, if
    !         itnlim > 0, min(options%localSize,options%itnlim,m,n).
    ! See also options%orthog_tol

    integer(ip) :: stop_test = 1 ! set up data required for
      !              using the Papez and Tichy stopping test.
      ! The user performs the test for convergence

    ! stop_test = 2: this is as for stop_test = 1 except that
      ! the code performs the Papez and Tichy stopping test.
      ! For this, the user is required to supply an estimate of ||A||_2 using
      ! the optional input argument anorm. If this is not present, stop_test=1
      ! is used. Note: stop_test = 2 should NOT be used if precon = 0
      ! and split preconditioning is used (so that A is replaced by A*L^{-1}).
      ! In this case, the user must test for convergence (otherwise, the
      ! stopping test is for the preconditioned problem).

    ! if stop_test = 0, then this data is not set up (saves memory and work but
      ! user then has to decide how to terminate the computation).

    ! other values treated as stop_test = 1

    real(wp) :: delta = sqrt(epsilon(1.0_wp)) ! convergence tolerance 
      ! if stop_test = 2

    integer(ip) :: orthog_choice = 0 ! controls reorthog strategy
    ! 0 = no reorthogonalisation 
    ! 1 = reorthogonalise using MGS against previous options%localSize vectors
    ! 2 = reorthogonalise using MGS against first options%localSize vectors
    !-1 = reorthogonalise using MGS2 against previous options%localSize vectors
    !-2 = reorthogonalise using MGS2 against first options%localSize vectors
    ! All other values are treated as zero
    ! MGS  = modified Gram Schmidt
    ! MGS2 = two applications of MGS (can help but can add significant overhead)

    real(wp) :: orthog_tol = zero ! if orthog_tol > zero then
      ! selective reorthogonalization (ie only reorthog v against a
      ! previous vector if inner product between them is greater
      ! than orthog_tol). MGS2 is not used with this option (so abs(orthog_choice)
      ! is used).

    real(wp) :: stagnation = zero ! stagnation occurs
    ! if the update to x on some iteration is less than options%stagnation

    ! The following are used by Papez and Tichy stopping (these values were used
    ! in their paper). Only accessed on first call. Values .le. zero are reset
    ! to the default.
    real(wp) :: tau = 0.25_wp   ! 
    real(wp) :: tol = 0.0001_wp ! 

  end type lsqr_options

! ------------------------------------------------------------------

  ! information data type
  type :: lsqr_inform

    integer(ip) :: flag ! Gives reason for termination 
    ! negative indicates failure

    ! lsqr_stop_m_oor        :  n < 1 or m < 1
    ! lsqr_stop_allocation   :  An array allocation failed.
    ! lsqr_stop_deallocation :  An array deallocation failed.
    ! lsqr_stop_itnlim       :  iteration limit reached
    !                           x holds the current solution.

    ! warning flags(explained below)
    ! lsqr_stop_anorm     : options%stop_test = 2 only.
    !    Either the user did not supply anorm or anorm.le.0
    ! lsqr_stop_x0        : x = 0  is the exact solution.
    !    No iterations were performed.
    ! lsqr_stop_stagnate  : computation of LS solution has stagnated.
    !    User should test for convergence (options%stop_test.ne.2)
    !    and then terminate (even if convergence test not satisfied)

    integer(ip) :: itn    ! The number of iterations performed

    integer :: stat   ! Fortran stat parameter

    real(wp) :: estim ! set if options%stop_test = 1 or 2. In this case, holds 
    ! backward error estimate that is required by Papez-Tichy stopping test.

    real :: time_RO ! accumulate the time for performing reorthogonalization.

  end type lsqr_inform

! ------------------------------------------------------------------

!  Derived type to ensure local variables are saved safely in the
!  reverse communication
  type :: lsqr_keep
    private

    logical     :: damped, localVQueueFull

    integer(ip) :: branch = 0
    integer(ip) :: flag = 0
    integer(ip) :: itn_count ! iteration count
    integer(ip) :: itnlim
    integer(ip) :: localPointer, localVecs, orthog_choice
    integer(ip) :: precon
    integer(ip) :: stop_test  ! copy of options%stop_test
    integer(ip) :: test_frequency ! testing frequency

    real(wp)    :: alpha, beta, rho, rhobar
    real(wp)    :: backward
    real(wp)    :: anorm, bnorm, xnorm
    real(wp)    :: orthog_tol ! copy of options%orthog_tol
    real(wp)    :: damp
    real(wp)    :: phibar
    real(wp)    :: stagnation ! copy of options%stagnation
    real(wp)    :: t1, t2

    real(wp), allocatable :: w(:)
    real(wp), allocatable :: localV(:,:)

    ! These are only needed if stop_test = 1 or 2
    !  (Papez and Tichy stopping condition)
    integer(ip) :: d, ell, estimend
    real(wp)    :: tau, tol ! copies of options%tau, options%tol

    real(wp), allocatable, dimension(:) :: delay
    real(wp), allocatable, dimension(:) :: delta
    real(wp), allocatable, dimension(:) :: curve
    real(wp), allocatable, dimension(:) :: estim

  end type lsqr_keep

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lsqr_double (precon, action, m, n, u, v, z, x, keep, &
        options, inform, damp, anorm)

    integer(ip), intent(in)  :: precon  ! This parameter controls
    ! the use of preconditioning. It must be set on the first call
    ! to one of the following values:

    ! 0:  No preconditioning. Also can be used for right preconditioning by
    !     replacing products with A by products with A* P^{-T} 
    !     where the preconditioner is M = PP^T. In this case,
    !     the least squares solution must be recovered by the user
    !     as x_LS = P^T x, where x is returned by lsqr (if testing is
    !     to be done for the unpreconditioned problem, this will need to be
    !     performed each time action = 4 is returned). 
    !     The advantage compared to precon = 1 is that the array z is not
    !     needed (and so can be of size 1).
    ! 1:  Right (split) preconditioning. Requires the preconditioner
    !     is in factored form M = PP^T eg if A^T*A = LL^T (Cholesky fact)
    !     can use P = L
    ! 2:  Left preconditioning of normal equations. M does not need to
    !     be in factored form.

    ! Note: precon is only accessed on the first call (the call with 
    ! action = 0). If the user wishes to try different options, then
    ! the computation must be restarted.

    integer(ip), intent(inout) :: action ! This parameter controls
    ! the action. On initial entry, action must be set to 0
    ! and the initial guess for the solution must be in the array x
    ! (if no guess is available, the user should set x(1:n) = 0).
    ! On subsequent calls, action should be unchanged by the user.
    ! On return, action determines what the user is required to do.
    ! Possible values and their consequences are:
    ! 
    ! 0. Computation has terminated.
    !    action = 0 is returned if an error has been encountered.
    !    Also returned if options%stop_test = 2 and convergence 
    !    has been achieved
    !
    ! 1. precon = 0 or 2 : the user must compute v = v + A^T *u
    !    precon = 1 : the user must compute v = v + P^{-1} A^T *u 
    ! The user must then re-call the subroutine without any other arguments. 
    ! The vectors u and v are available in the arrays u and v (see below).
    ! This action is returned once per iteration.
    
    ! 2. precon = 0 : the user must compute u = u + A *v
    !    precon = 1 or 2 : the user must compute u = u + A *z 
    ! The user must then re-call the subroutine without any other arguments.
    ! The vectors u, v and z are available in the arrays u, v and z (see below).
    ! This action is returned once per iteration.

    ! 3. precon = 1 : the user must compute z = P^{-T}*v
    !    precon = 2 : the user must compute z = M^{-1}*v
    ! The user must then re-call the subroutine without any other arguments.
    ! The vectors z and v are available in the arrays z and v (see below). 
    ! This action is returned once per iteration.

    ! 4:  The user may test for convergence.
    !     Re-call the subroutine WITHOUT alterning any of the arguments.
    !     This return does not happen if options%stop_test = 2

    integer(ip), intent(in)  :: m  ! the number of rows in A.

    integer(ip), intent(in)  :: n  ! the number of columns in A.

    real(wp),    intent(inout) :: u(m) ! Prior to the first call, u
    ! must hold the rhs vector b. u is then the vector used to communicate u 
    ! when further action is required (see action, above).

    real(wp),    intent(inout) :: v(n) ! The vector used to communicate v 
    !  when further action is required (see action, above).

    real(wp),    intent(inout) :: z(:) ! The vector used to communicate z
    !  when further action is required (see action, above). 
    !  z is not accessed if precon = 0 (so can be size 1 in this case;
    !  otherwise it should be length n)

    real(wp),    intent(inout) :: x(n) ! Does not need to be set by the user.
    ! Returns the computed solution x

    type ( lsqr_keep ), intent( inout) :: keep 
    ! A variable of type lsqr_keep that is used to
    ! preserve internal data between reverse-communication
    ! calls. The components are private.

    type ( lsqr_options ), intent( inout) :: options
    ! A variable of type lsqr_options that is used to control the options.
    ! It should not be altered by the user between calls.

    type ( lsqr_inform ), intent( inout) :: inform
    ! A variable of type lsqr_info_type that is used to hold information.
    ! It should not be altered by the user (so as to preserve the
    ! information between calls).

    real(wp),intent(in), optional  :: damp ! The damping parameter.
    ! Only accessed on the first call. damp should be > 0
    ! Note: The work per iteration and the storage needed
    ! by LSQR are the same for all values of damp.

    real(wp),intent(in), optional  :: anorm ! The norm of A.
    ! must be present and hold an estimate of the 2 norm of A
    ! if options%stop_test = 2


    real(wp) :: ddot, dnrm2
    ! Local arrays and variables
    integer     :: dst, st
    real(wp)    :: backward
    real(wp)    :: cs, cs1, phi, psi, rhbar1, rho, sn, sn1
    real(wp)    :: theta
    real(wp)    :: xnorm

    !-------------------------------------------------------------------

    ! on first call, initialize keep%branch and keep%flag
    ! And initial components of inform so they are not undefined.
    if (action.eq.0) then
       keep%branch = 0
       keep%flag   = 0

       inform%flag    = 0
       inform%itn     = 0
       inform%stat    = 0
       inform%estim   = -1.0_wp
       inform%time_RO   = 0.0
    end if

    ! Immediate return if we have already had an error 
    if (keep%flag.lt.0) then
       inform%flag = keep%flag
       action = 0
       return
    end if

    ! also terminate if we things have stagnated.
    if (keep%flag.ge.lsqr_stop_stagnate) then
       action = 0
       return
    end if

    ! on other calls, jump to the appropriate place after reverse-communication

    if (precon.eq.0) then
      ! LSQR without preconditioning (or split preconditioning with 
      ! A replaced by AP^{-T})
      select case ( keep%branch )
        case ( 1 )
          go to 20
        case ( 2 )
          go to 40
        case ( 3 )
          go to 50
        case ( 4 )
          go to 60
      end select

    else if (precon.eq.1) then
      ! modified RP-LSQR
      select case ( keep%branch )
        case ( 1 )
          go to 120
        case ( 2 )
          go to 140
        case ( 3 )
          go to 150
        case ( 4 )
          go to 160
        case ( 5 )
          go to 125
        case ( 6 )
          go to 170
      end select

    else if (precon.eq.2) then
      ! LP-LSQR
      select case ( keep%branch )
        case ( 1 )
          go to 220
        case ( 2 )
          go to 240
        case ( 3 )
          go to 250
        case ( 4 )
          go to 280
        case ( 5 )
          go to 225
        case ( 6 )
          go to 270
      end select
    end if

    ! Initialize.
    keep%damped  = .false.
    if (present(damp)) then
       keep%damp = damp
       keep%damped = keep%damp > zero
    end if

    keep%itnlim = options%itnlim
    if (options%itnlim < 1) keep%itnlim = n

    keep%stop_test = options%stop_test
    if (keep%stop_test.eq.2) then
      if (.not. present(anorm)) then
         keep%stop_test = 1
         keep%flag = lsqr_stop_anorm
      else if (anorm .le. zero) then
         keep%stop_test = 1
         keep%flag = lsqr_stop_anorm
      else
         keep%anorm = anorm
      end if
    else if (keep%stop_test .ne. 0) then
      keep%stop_test = 1
    end if
    inform%flag = keep%flag

    keep%itn_count = 0
    keep%test_frequency = options%test_frequency
    if (keep%test_frequency .lt. 1) keep%test_frequency = 1

    !!! orthogonalisation controls
    keep%orthog_choice = options%orthog_choice
    if (abs(keep%orthog_choice).gt.2) keep%orthog_choice = 0


    keep%localVecs = min(options%localSize,m,n)
    if (options%itnlim > 0) &
        keep%localVecs = min(keep%localVecs, options%itnlim)
    if (keep%orthog_choice.eq.0) keep%localVecs = 0

    keep%orthog_tol = max(zero, options%orthog_tol)
    ! Do not use MGS2 if we have orthog_tol non zero
    if (keep%orthog_tol.ne.zero) keep%orthog_choice = abs(keep%orthog_choice)
    !!!!

    keep%stagnation = max(zero, options%stagnation)

    keep%tau = options%tau
    if (keep%tau .le. zero) keep%tau = 0.25_wp

    keep%tol = options%tol
    if (keep%tol .le. zero) keep%tol = 0.0001_wp

    st  = 0 
    dst = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! quick check of m and n
    if (n.lt.1 .or. m.lt.1) then
       keep%flag   = lsqr_stop_m_oor
       inform%flag = lsqr_stop_m_oor
       return
    end if
 
    ! if precon is out-of-range, run without preconditioning
    keep%precon = precon
    if (precon.lt.1 .or. precon.gt.2) keep%precon = 0

    ! allocate workarrays

    if (.not. allocated(keep%w)) then
      allocate( keep%w(n), stat = st)
      if (st /= 0) go to 10
    else if (size(keep%w).lt.n) then
      deallocate (keep%w, stat = dst)
      if (dst /= 0) then
        go to 10
      else
        allocate( keep%w(n), stat = st)
        if (st /= 0) go to 10
      end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set up buffer for optional reorthogonalization
    if (keep%localVecs > 0) then
      if (.not. allocated(keep%localV)) then
        allocate( keep%localV(n,keep%localVecs), stat = st)
        if (st /= 0) go to 10
      else if (size(keep%localV,1).lt.n .or. &
        size(keep%localV,2).lt.keep%localVecs) then
        deallocate (keep%localV, stat = dst)
        if (dst /= 0) go to 10
        allocate( keep%localV(n,keep%localVecs), stat = st)
        if (st /= 0) go to 10
      end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set up for Papez and Tichy stopping
    if (keep%stop_test.ne.0) then
      keep%ell = 1
      keep%d = 0
      keep%estimend = -1

      if (allocated(keep%curve)) deallocate(keep%curve, stat = dst)
      if (dst /= 0) go to 10
      if (allocated(keep%delay)) deallocate(keep%delay, stat = dst)
      if (dst /= 0) go to 10
      if (allocated(keep%delta)) deallocate(keep%delta, stat = dst)
      if (dst /= 0) go to 10
      if (allocated(keep%estim)) deallocate(keep%estim, stat = dst)
      if (dst /= 0) go to 10
 
      allocate (keep%curve(keep%itnlim), keep%delay(keep%itnlim),  &
           keep%delta(keep%itnlim), keep%estim(keep%itnlim), stat = st)
      if (st /= 0) go to 10
    end if

10 continue
   if (st.ne.0) then
      inform%stat = st
      inform%flag = lsqr_stop_allocation
      return
   else if (dst.ne.0) then
      inform%stat = dst
      inform%flag = lsqr_stop_deallocation
      return
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

    keep%bnorm   = dnrm2 (m, u, 1)
    if (keep%bnorm .eq. zero) then
      ! Exit if b = 0. Solution is x_{LS} = 0.
       x(1:n) = zero
       keep%flag   = keep%flag + lsqr_stop_x0
       inform%flag = keep%flag
       return
    end if

    !-------------------------------------------------------------------
    ! Set up the first vectors u_1 and v_1 for the bidiagonalization.
    ! These satisfy  beta*u = r_0 = b,  alpha*v = A^T *u.
    ! u holds r_0 = b - A*x_0. 
    !-------------------------------------------------------------------

    keep%beta   = dnrm2 (m, u, 1)

    call PREC(scal)(m, (one/keep%beta), u, 1)  ! u_1 = u / beta

    ! User must compute v = A^T * u (precon = 0 or 2),
    ! or v = P^{-1} A^T *u (precon = 1)

    ! initialise v = 0 (as user may have routine for computing v = v + A^T * u)
    v(1:n) = zero

    action = 1              
    keep%branch = 1
    return

    !===================================================================
    !===================================================================
    ! precon = 0 : This is the no preconditioning variant
    !              or split conditioning with A replaced by AP^{-T} 

 20 continue

    ! v = A^T *u_1 has been performed by the user

    keep%alpha = dnrm2 (n, v, 1)

    if (keep%alpha .eq. zero) then
       ! terminate A'b = 0
       action = 0
       return
    end if

    call PREC(scal)(n, (one/keep%alpha), v, 1)

    ! w_1 = v_1
    keep%w = v

    ! Initialization for local reorthogonalization.
    if (keep%localVecs > 0) then
       keep%localPointer    = 1
       keep%localVQueueFull = .false.
       keep%localV(1:n,1)   = v(1:n)
    end if

    ! Initialize variables for 1st iteration.
    keep%rhobar   = keep%alpha
    keep%phibar   = keep%beta

    !===================================================================
    ! Main iteration loop.
    !===================================================================

 30 continue

       keep%itn_count = keep%itn_count + 1
       inform%itn = keep%itn_count

       !----------------------------------------------------------------
       ! Perform the next step of the bidiagonalization to obtain the
       ! next beta, u, alpha, v.  These satisfy
       !     beta*u =    A*v  - alpha*u,
       !    alpha*v = (A)'*u  -  beta*v.
       !----------------------------------------------------------------
       call PREC(scal)(m,(- keep%alpha), u, 1)  ! - alpha*u_j

       ! Compute u = u + A*v
       action = 2
       keep%branch = 2
       return

 40    continue

       ! user has computed u_{j+1} = A*v_j - alpha *u_j. Normalise.
       keep%beta   = dnrm2 (m, u, 1)

       if (keep%beta > zero) then
          call PREC(scal)(m, (one/keep%beta), u, 1) ! u = u / beta
          if (keep%localVecs > 0) then
             ! Check whether to store v_j for local reorthog'n of v_{j+1}
             if (abs(keep%orthog_choice).eq.1 .or.       &
                keep%itn_count .lt. keep%localVecs)      & 
                call localVEnqueue(n, v, keep)
          end if  

          call PREC(scal)(n, (-keep%beta), v, 1) ! v = -beta *v_j

          ! Compute v = v + A'*u
          action = 1
          keep%branch = 3
          return
       end if

 50    continue

       ! user has computed v_{j+1} = A^T*u_{j+1} - beta *v_j.
       if (keep%beta > zero) then
          if (keep%localVecs > 0) & ! Local reorthogonalization of new v.
             call localVOrtho(n, v, keep, inform%time_RO)

          keep%alpha  = dnrm2 (n, v, 1)
          if (keep%alpha > zero) call PREC(scal)(n, (one/keep%alpha), v, 1)

       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the damping parameter.
       ! This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
       !----------------------------------------------------------------
       rhbar1 = keep%rhobar
       psi    = zero
       if (keep%damped) then
          rhbar1 = d2norm(keep%rhobar, keep%damp)
          cs1    = keep%rhobar/rhbar1
          sn1    = keep%damp  /rhbar1
          psi    = sn1 * keep%phibar
          keep%phibar = cs1 * keep%phibar
       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the subdiagonal element (beta)
       ! of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
       !----------------------------------------------------------------
       rho    =   d2norm(rhbar1, keep%beta)
       cs     =   rhbar1 /rho
       sn     =   keep%beta /rho
       theta  =   sn*keep%alpha
       keep%rhobar = - cs*keep%alpha
       phi         =   cs*keep%phibar
       keep%phibar = sn*keep%phibar

       !----------------------------------------------------------------
       ! Update  x
       ! ---------------------------------------------------------------
       keep%t1     = phi/rho
       keep%t2     = - theta/rho
       x(1:n)      = keep%t1* keep%w(1:n) + x(1:n)

       ! the following enables Papez and Tichy stopping test
       if (keep%stop_test.ne.0) then
          call adaptive(keep%itn_count,phi*phi,keep)
          if (keep%estimend > 0) then
             inform%estim = keep%estim(keep%estimend)
          end if
       end if

       !----------------------------------------------------------------
       ! Test for convergence. And also check for stagnation.
       !----------------------------------------------------------------

       if (keep%itn_count >= keep%itnlim) then ! iteration count exceeded
          inform%flag = lsqr_stop_itnlim
          keep%flag   = lsqr_stop_itnlim
          action = 0
          return
       end if

       if (abs(keep%t1).le.keep%stagnation) then
          ! Stagnated.
          if (keep%stop_test.ne.2) then
            ! The user should check for convergence.
            if (keep%flag .le. lsqr_stop_stagnate) then
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
            end if
            action = 4
            keep%branch = 4
            return
          else
            ! Perform Papez and Tichy testing
            xnorm = dnrm2 (n, x, 1)
            backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)

            if (backward > 0 .and. backward < options%delta) then
              ! convergence achieved
              action = 0
              return
            else
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
              action = 0
              return
            end if
          end if
       end if

       ! see if it is time to return to test convergence
       if (mod(keep%itn_count,keep%test_frequency) == 0) then
          if (keep%stop_test.ne.2) then
            action = 4
            keep%branch = 4
            return
          end if
       end if
       if (keep%stop_test.eq.2) then
          ! Perform Papez and Tichy testing (do this every iteration
          ! as involves little extra work)
          xnorm = dnrm2 (n, x, 1)
          backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)
          if (backward > 0 .and. backward < options%delta) then
            ! convergence achieved
            action = 0
            return
          end if
       end if

  60   continue

       !----------------------------------------------------------------
       ! Update  w
       ! ---------------------------------------------------------------

       keep%w(1:n) = keep%t2* keep%w(1:n) + v(1:n)

      ! cycle for next iteration       
       go to 30

    !===================================================================
    !===================================================================
    ! This is the modified RP-LSQR preconditioning variant (precon = 1)

120 continue

    ! v =  P^{-1} A^T *u has been performed by the user

    keep%alpha = dnrm2 (n, v, 1)

    if (keep%alpha .eq. zero) then
       ! terminate P^{-1} A^T b = 0
       action = 0
       return
    end if

    call PREC(scal)(n, (one/keep%alpha), v, 1)

    !----------------------------------------------------------------
    ! user to compute z_1 = P^{-T} v_1
    !----------------------------------------------------------------

    action = 3
    keep%branch = 5
    return

 125 continue

    ! d_1 = z_1
    keep%w = z

    ! Initialization for local reorthogonalization.
    if (keep%localVecs > 0) then
       keep%localPointer    = 1
       keep%localVQueueFull = .false.
       keep%localV(1:n,1)   = v(1:n)
    end if

    ! Initialize variables for 1st iteration.
    keep%rhobar   = keep%alpha
    keep%phibar   = keep%beta

    !===================================================================
    ! Main iteration loop.
    !===================================================================

130 continue

       keep%itn_count = keep%itn_count + 1
       inform%itn = keep%itn_count

       !----------------------------------------------------------------
       ! Perform the next step of the bidiagonalization to obtain the
       ! next beta, u, alpha, v.  These satisfy
       !     beta*u =    A*z  - alpha*u,
       !    alpha*v = P^{-1} A^T*u  -  beta*v.
       !----------------------------------------------------------------
       call PREC(scal)(m,(-keep%alpha), u, 1)  ! -alpha*u_j

       ! Compute u = u + A*z
       action = 2
       keep%branch = 2
       return

140    continue

       ! user has computed u_{j+1} = A*z_j - alpha *u_j. Normalise.
       keep%beta   = dnrm2 (m, u, 1)

       if (keep%beta > zero) then
          call PREC(scal)(m, (one/keep%beta), u, 1) ! u = u / beta
          if (keep%localVecs > 0) then
             ! Check whether to store v_j for local reorthog'n of v_{j+1}
             if (abs(keep%orthog_choice).eq.1 .or.     &
                keep%itn_count .lt. keep%localVecs)    & 
                call localVEnqueue(n, v, keep)
          end if 

          call PREC(scal)(n, (-keep%beta), v, 1) ! v = -beta *v_j

          ! Compute v = v + P^{-1} A^T*u
          action = 1
          keep%branch = 3
          return
       end if

150    continue

       ! user has computed v_{j+1} = A^T*u_{j+1} - beta *v_j.
       if (keep%beta > zero) then

          if (keep%localVecs > 0) & ! Local reorthogonalization of new v.
             call localVOrtho(n, v, keep, inform%time_RO)

          keep%alpha  = dnrm2 (n, v, 1)
          if (keep%alpha > zero) call PREC(scal)(n, (one/keep%alpha), v, 1)

       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the damping parameter.
       ! This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
       !----------------------------------------------------------------
       rhbar1 = keep%rhobar
       psi    = zero
       if (keep%damped) then
          rhbar1 = d2norm(keep%rhobar, keep%damp)
          cs1    = keep%rhobar/rhbar1
          sn1    = keep%damp  /rhbar1
          psi    = sn1 * keep%phibar
          keep%phibar = cs1 * keep%phibar
       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the subdiagonal element (beta)
       ! of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
       !----------------------------------------------------------------
       rho    =   d2norm(rhbar1, keep%beta)
       cs     =   rhbar1 /rho
       sn     =   keep%beta /rho
       theta  =   sn*keep%alpha
       keep%rhobar = - cs*keep%alpha
       phi         =   cs*keep%phibar
       keep%phibar = sn*keep%phibar

       !----------------------------------------------------------------
       ! Update  x
       ! ---------------------------------------------------------------
       keep%t1  =     phi/rho
       keep%t2  = - theta/rho

       x(1:n)   = keep%t1* keep%w(1:n) + x(1:n)

       ! the following enables Papez and Tichy stopping test
       if (keep%stop_test.ne.0) then
          call adaptive(keep%itn_count,phi*phi,keep)
          if (keep%estimend > 0) & 
             inform%estim = keep%estim(keep%estimend)   
       end if

       !----------------------------------------------------------------
       ! Test for convergence. And also check for stagnation.
       !----------------------------------------------------------------

       if (keep%itn_count >= keep%itnlim) then ! iteration count exceeded
          inform%flag = lsqr_stop_itnlim
          keep%flag   = lsqr_stop_itnlim
          action = 0
          return
       end if

       if (abs(keep%t1).le.keep%stagnation) then
          ! Stagnated.
          if (keep%stop_test.ne.2) then
            ! The user should check for convergence.
            if (keep%flag .le. lsqr_stop_stagnate) then
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
            end if
            action = 4
            keep%branch = 4
            return
          else
            ! Perform Papez and Tichy testing
            xnorm = dnrm2 (n, x, 1)
            backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)
            if (backward > 0 .and. backward < options%delta) then
              ! convergence achieved
              action = 0
              return
            else
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
              action = 0
              return
            end if
          end if
       end if

       ! see if it is time to return to test convergence
       if (mod(keep%itn_count,keep%test_frequency) == 0) then
          if (keep%stop_test.ne.2) then
            action = 4
            keep%branch = 4
            return
          end if
       end if
       if (keep%stop_test.eq.2) then
          ! Perform Papez and Tichy testing (do this every iteration
          ! as involves little extra work)
          xnorm = dnrm2 (n, x, 1)
          backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)
          if (backward > 0 .and. backward < options%delta) then
            ! convergence achieved
            action = 0
            return
          end if
       end if

 160   continue
       !----------------------------------------------------------------
       ! User to compute z = P^{-T} *v
       !----------------------------------------------------------------
       action = 3
       keep%branch = 6
       return

 170   continue

       !----------------------------------------------------------------
       ! Update d
       ! ---------------------------------------------------------------

       keep%w(1:n) = keep%t2* keep%w(1:n) + z(1:n)

      ! cycle for next iteration       
       go to 130


    !===================================================================
    !===================================================================
    ! This is the LP-LSQR (left) preconditioning variant (precon = 2). 
    ! Uses M^{-1} (factored form not needed

220 continue

    ! v =  A^T *u has been performed by the user

    !----------------------------------------------------------------
    ! user to compute z_1 = M^{-1} v_1
    !----------------------------------------------------------------

    action = 3
    keep%branch = 5
    return

 225 continue

    keep%alpha = sqrt( ddot (n, v, 1, v, 1))

    if (keep%alpha .eq. zero) then
       ! terminate A^T b = 0
       action = 0
       return
    end if

    call PREC(scal)(n, (one/keep%alpha), v, 1)
    call PREC(scal)(n, (one/keep%alpha), z, 1)

    ! d_1 = z_1
    keep%w = z

    ! Initialization for local reorthogonalization.
    if (keep%localVecs > 0) then
       keep%localPointer    = 1
       keep%localVQueueFull = .false.
       keep%localV(1:n,1)   = v(1:n)
    end if

    ! Initialize variables for 1st iteration.
    keep%rhobar   = keep%alpha
    keep%phibar   = keep%beta

    !===================================================================
    ! Main iteration loop.
    !===================================================================

230 continue

       keep%itn_count = keep%itn_count + 1
       inform%itn = keep%itn_count

       !----------------------------------------------------------------
       ! Perform the next step of the bidiagonalization to obtain the
       ! next beta, u, alpha, v.  These satisfy
       !     beta*u =    A*z  - alpha*u,
       !    alpha*v = (A)'*u  -  beta*v.
       !----------------------------------------------------------------
       call PREC(scal)(m,(- keep%alpha), u, 1)  ! - alpha*u_j
       ! Compute u = u + A*z
       action = 2
       keep%branch = 2
       return

240    continue

       ! user has computed u_{j+1} = A*z_j - alpha *u_j. Normalise.
       keep%beta   = dnrm2 (m, u, 1)

       if (keep%beta > zero) then
          call PREC(scal)(m, (one/keep%beta), u, 1) ! u = u / beta
          if (keep%localVecs > 0) then
             ! Check whether to store v_j for local reorthog'n of v_{j+1}
             if (abs(keep%orthog_choice).eq.1 .or.    &
                keep%itn_count .lt. keep%localVecs)   & 
                call localVEnqueue(n, v, keep)
          end if 

          call PREC(scal)(n, (-keep%beta), v, 1) ! v = -beta *v_j

          ! Compute v = v + (A)'*u
          action = 1
          keep%branch = 3
          return
       end if

250    continue

       !----------------------------------------------------------------
       ! User to compute z = M^{-1} *v
       !----------------------------------------------------------------
       if (keep%localVecs > 0) & ! Local-reorthogonalization of new v.
          call localVOrtho(n, v, keep, inform%time_RO)

       action = 3
       keep%branch = 6
       return

 270 continue

       ! user has computed v_{j+1} = A^T*u_{j+1} - beta *v_j and 
       ! z_{j+1} = M^{-1} *v_{j+1}
       if (keep%beta > zero) then

          keep%alpha  = sqrt( ddot(n, v, 1, z, 1))
          if (keep%alpha > zero) then
            call PREC(scal)(n, (one/keep%alpha), v, 1)
            call PREC(scal)(n, (one/keep%alpha), z, 1)
          end if
       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the damping parameter.
       ! This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
       !----------------------------------------------------------------
       rhbar1 = keep%rhobar
       psi    = zero
       if (keep%damped) then
          rhbar1 = d2norm(keep%rhobar, keep%damp)
          cs1    = keep%rhobar/rhbar1
          sn1    = keep%damp  /rhbar1
          psi    = sn1 * keep%phibar
          keep%phibar = cs1 * keep%phibar
       end if

       !----------------------------------------------------------------
       ! Use a plane rotation to eliminate the subdiagonal element (beta)
       ! of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
       !----------------------------------------------------------------
       rho    =   d2norm(rhbar1, keep%beta)
       cs     =   rhbar1 /rho
       sn     =   keep%beta /rho
       theta  =   sn*keep%alpha
       keep%rhobar = - cs*keep%alpha
       phi         =   cs*keep%phibar
       keep%phibar = sn*keep%phibar

       !----------------------------------------------------------------
       ! Update  x
       ! ---------------------------------------------------------------
       keep%t1  =     phi/rho
       keep%t2  = - theta/rho

       x(1:n)   = keep%t1* keep%w(1:n) + x(1:n)

       ! the following enables Papez and Tichy stopping test
       if (keep%stop_test.ne.0) then
          call adaptive(keep%itn_count,phi*phi,keep)
          if (keep%estimend > 0) & 
             inform%estim = keep%estim(keep%estimend)   
       end if

       !----------------------------------------------------------------
       ! Test for convergence. And also check for stagnation.
       !----------------------------------------------------------------

       if (keep%itn_count >= keep%itnlim) then ! iteration count exceeded
          inform%flag = lsqr_stop_itnlim
          keep%flag   = lsqr_stop_itnlim
          action = 0
          return
       end if

       if (abs(keep%t1).le.keep%stagnation) then
          ! Stagnated.
          if (keep%stop_test.ne.2) then
            ! The user should check for convergence.
            if (keep%flag .le. lsqr_stop_stagnate) then
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
            end if
            action = 4
            keep%branch = 4
            return
          else
            ! Perform Papez and Tichy testing
            xnorm = dnrm2 (n, x, 1)
            backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)
            if (backward > 0 .and. backward < options%delta) then
              ! convergence achieved
              action = 0
              return
            else
              keep%flag   = keep%flag + lsqr_stop_stagnate
              inform%flag = keep%flag
              action = 0
              return
            end if
          end if
       end if

       ! see if it is time to return to test convergence
       if (mod(keep%itn_count,keep%test_frequency) == 0) then
          if (keep%stop_test.ne.2) then
            action = 4
            keep%branch = 4
            return
          end if
       end if
       if (keep%stop_test.eq.2) then
          ! Perform Papez and Tichy testing (do this every iteration
          ! as involves little extra work)
          xnorm = dnrm2 (n, x, 1)
          backward = inform%estim /(keep%anorm*xnorm + keep%bnorm)
      ! write (*,*) xnorm, backward
          if (backward > 0 .and. backward < options%delta) then
            ! convergence achieved
            action = 0
            return
          end if
       end if

 280   continue
       !----------------------------------------------------------------
       ! Update d
       ! ---------------------------------------------------------------

       keep%w(1:n) = keep%t2* keep%w(1:n) + z(1:n)

      ! cycle for next iteration       
       go to 230

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function d2norm( a, b )

      real(wp)             :: d2norm
      real(wp), intent(in) :: a, b

      !-------------------------------------------------------------------
      ! d2norm returns sqrt( a**2 + b**2 )
      ! with precautions to avoid overflow.
      !-------------------------------------------------------------------

      intrinsic            :: abs, sqrt
      real(wp)             :: scale

      d2norm = zero
      scale = abs(a) + abs(b)
      if (scale > zero) d2norm = scale*sqrt((a/scale)**2 + (b/scale)**2)

    end function d2norm

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine localVEnqueue(n, v, keep)

    integer(ip), intent(in) :: n
    real(wp), intent(in)    :: v(n)
    type ( lsqr_keep ), intent(inout) :: keep

      ! Store v into the circular buffer keep%localV.

      if (keep%localPointer < keep%localVecs) then
         keep%localPointer = keep%localPointer + 1
      else
         keep%localPointer = 1
         keep%localVQueueFull = .true.
      end if
      keep%localV(1:n,keep%localPointer) = v

    end subroutine localVEnqueue

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine localVOrtho(n, v, keep, time)

   ! Perform local reorthogonalization of current v
   ! against previously computed v's (which are orthonormal
   ! so that ddot(n,keep%localV(1:n,j),1,keep%localV(1:n,j),1) = 1
   ! and does not need computing).
   ! Modified Gram Schmidt reorthog. is used (MGS2 if keep%orthog_choice < 0)
   ! Notes: Found it is more efficient to use blas ddot (rather than dot_product).
   ! Using precompiled blas can lead to different performance (LSQR iteration counts)

    integer(ip), intent(in) :: n
    real(wp), intent(inout) :: v(n)
    type ( lsqr_keep ), intent(in) :: keep
    real, intent(inout) :: time

    real(wp)    :: d, ddot
    integer(ip) :: j, limit
    integer(long) :: time1, time2, rate_t

      ! set limit to be number of vectors we are orthogonalising against
      if (keep%localVQueueFull) then
         limit = keep%localVecs
      else
         limit = keep%localPointer
      end if

      ! time the cost of RO
      call system_clock(time1, rate_t)

      if (keep%orthog_tol > zero) then
        ! selective reorthogonalisation
        do j = 1, limit
          d = ddot(n,v,1,keep%localV(1:n,j),1)
          if (d .gt. keep%orthog_tol) call PREC(axpy)(n,-d,keep%localV(1:n,j),1,v,1)
        end do
      else
        do j = 1, limit
          d = ddot(n,v,1,keep%localV(1:n,j),1)
          call PREC(axpy)(n,-d,keep%localV(1:n,j),1,v,1)
          if (keep%orthog_choice .lt. 0) then
            ! second application of MGS
            d = ddot(n,v,1,keep%localV(1:n,j),1)
            call PREC(axpy)(n,-d,keep%localV(1:n,j),1,v,1)
          end if
        end do
      end if

      call system_clock(time2)

      time = time + (time2 - time1)/real(rate_t)

    end subroutine localVOrtho
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine lsqr_double

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lsqr_free_double (keep, stat)

  ! Routine to deallocate components of keep. Problem is indicated
  ! by nonzero stat value. No printing.

    type ( lsqr_keep ), intent( inout) :: keep
    integer(ip), intent(out) :: stat ! set to Fortran stat parameter
                                     ! Non zero only if a deallocation fails
    integer(ip) :: st

    stat = 0

    deallocate (keep%w, stat=st)
    if (st.ne.0) stat = st

    if (allocated(keep%localV)) then
       deallocate (keep%localV, stat=st)
       if (st.ne.0) stat = st
    end if

    ! If Papez and Tichy stopping was used, additional arrays to be dellocated
    if (keep%stop_test.ne.0) then
       deallocate (keep%delay, stat=st)
       if (st.ne.0) stat = st
       deallocate (keep%delta, stat=st)
       if (st.ne.0) stat = st
       deallocate (keep%curve, stat=st)
       if (st.ne.0) stat = st
       deallocate (keep%estim, stat=st)
       if (st.ne.0) stat = st
    end if

    keep%flag = 0 ! make sure error flag is reset as it is tested on the first
                  ! call to lsqr and we don't want to pick up an old value
    
  end subroutine lsqr_free_double

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine adaptive(k,add,keep)

      ! This is to get the information for the Papez and Tichy stopping
      ! criteria

      integer,  intent(in) :: k ! iteration number
      real(wp), intent(in) :: add

      type(lsqr_keep), intent(inout) :: keep
           ! components that are changed:
           ! curve(1:k)
           ! d
           ! delay(k)
           ! ell
           ! estim
           ! estimend

      integer  :: ell 
      real(wp) :: den, num, s

      ell = keep%ell

      keep%curve(1:k-1) = keep%curve(1:k-1) + add
      keep%curve(k) = add

      if (k > 1) then

        call find_s(k, keep%curve, keep%delta, ell, &
            keep%tol, s)

        num = s * add
        den = sum(keep%delta(ell:k-1))
              
        do while (keep%d >= 0 .and. num <= keep%tau*den)
          keep%delay(ell) = keep%d
          keep%estim(ell) = sqrt(den + add)
          keep%estimend = ell
          ell = ell + 1 
          keep%d = keep%d - 1
          den = sum(keep%delta(ell:k-1))
        end do
        ! update for next call
        keep%d = keep%d + 1
        keep%ell = ell
      end if
      keep%delta(k) = add

      end subroutine adaptive

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine find_s(k,curve,delta,ell,tol,s)

!      This is used by subroutine adaptive

      integer,  intent(in)  :: k,ell
      real(wp), intent(in)  :: curve(k),delta(k)
      real(wp), intent(in)  :: tol
      real(wp), intent(out) :: s
!
      integer  :: i,ind
      real(wp) :: temp
!      
      ind = 1
      do i = k,1,-1
        temp = curve(ell)/curve(i)
        if (temp <= tol) ind = i
      end do      

      s = zero  
      do i = ind, k-1
        s = max(s, curve(i)/delta(i))  
      end do    
!      
      end subroutine find_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lsqr_reverse_double
