! nlls_module :: a nonlinear least squares solver

module nlls_module

  implicit none

  integer, parameter :: wp = kind(1.0d0)
  integer, parameter :: long = selected_int_kind(8)
  real (kind = wp), parameter :: tenm5 = 1.0e-5
  real (kind = wp), parameter :: tenm8 = 1.0e-8
  real (kind = wp), parameter :: epsmch = epsilon(1.0_wp)
  real (kind = wp), parameter :: hundred = 100.0
  real (kind = wp), parameter :: ten = 10.0
  real (kind = wp), parameter :: point9 = 0.9
  real (kind = wp), parameter :: zero = 0.0
  real (kind = wp), parameter :: one = 1.0
  real (kind = wp), parameter :: two = 2.0
  real (kind = wp), parameter :: half = 0.5
  real (kind = wp), parameter :: sixteenth = 0.0625

  
  TYPE, PUBLIC :: NLLS_control_type
     
!   error and warning diagnostics occur on stream error 
     
     INTEGER :: error = 6

!   general output occurs on stream out

     INTEGER :: out = 6

!   the level of output required. <= 0 gives no output, = 1 gives a one-line
!    summary for every iteration, = 2 gives a summary of the inner iteration
!    for each iteration, >= 3 gives increasingly verbose (debugging) output

     INTEGER :: print_level = 0

!   any printing will start on this iteration

!$$     INTEGER :: start_print = - 1

!   any printing will stop on this iteration

!$$     INTEGER :: stop_print = - 1

!   the number of iterations between printing

!$$     INTEGER :: print_gap = 1

!   the maximum number of iterations performed

     INTEGER :: maxit = 100

!   removal of the file alive_file from unit alive_unit terminates execution

!$$     INTEGER :: alive_unit = 40
!$$     CHARACTER ( LEN = 30 ) :: alive_file = 'ALIVE.d'

!   non-monotone <= 0 monotone strategy used, anything else non-monotone
!     strategy with this history length used

!$$     INTEGER :: non_monotone = 1

!   specify the model used. Possible values are
!
!      0  dynamic (*not yet implemented*)
!      1  first-order (no Hessian)
!      2  second-order (exact Hessian)
!      3  barely second-order (identity Hessian)
!      4  secant second-order (sparsity-based)
!      5  secant second-order (limited-memory BFGS, with %lbfgs_vectors history)
!      6  secant second-order (limited-memory SR1, with %lbfgs_vectors history)

     INTEGER :: model = 1

!   specify the norm used. The norm is defined via ||v||^2 = v^T P v,
!    and will define the preconditioner used for iterative methods.
!    Possible values for P are
!
!     -3  user's own norm
!     -2  P = limited-memory BFGS matrix (with %lbfgs_vectors history)
!     -1  identity (= Euclidan two-norm)
!      0  automatic (*not yet implemented*)
!      1  diagonal, P = diag( max( Hessian, %min_diagonal ) )
!      2  banded, P = band( Hessian ) with semi-bandwidth %semi_bandwidth
!      3  re-ordered band, P=band(order(A)) with semi-bandwidth %semi_bandwidth
!      4  full factorization, P = Hessian, Schnabel-Eskow modification
!      5  full factorization, P = Hessian, GMPS modification (*not yet *)
!      6  incomplete factorization of Hessian, Lin-More'
!      7  incomplete factorization of Hessian, HSL_MI28
!      8  incomplete factorization of Hessian, Munskgaard (*not yet *)
!      9  expanding band of Hessian (*not yet implemented*)


     INTEGER :: nlls_method = 1

!   specify the method used to solve the trust-region sub problem
!      1 Powell's dogleg
!      ...


!$$     INTEGER :: norm = 1

!   specify the semi-bandwidth of the band matrix P if required

!$$     INTEGER :: semi_bandwidth = 5

!   number of vectors used by the L-BFGS matrix P if required

!$$     INTEGER :: lbfgs_vectors = 10

!   number of vectors used by the sparsity-based secant Hessian if required

!$$     INTEGER :: max_dxg = 100

!   number of vectors used by the Lin-More' incomplete factorization 
!    matrix P if required

!$$     INTEGER :: icfs_vectors = 10

!  the maximum number of fill entries within each column of the incomplete 
!  factor L computed by HSL_MI28. In general, increasing mi28_lsize improves
!  the quality of the preconditioner but increases the time to compute
!  and then apply the preconditioner. Values less than 0 are treated as 0

!$$     INTEGER :: mi28_lsize = 10

!  the maximum number of entries within each column of the strictly lower 
!  triangular matrix R used in the computation of the preconditioner by 
!  HSL_MI28.  Rank-1 arrays of size mi28_rsize *  n are allocated internally 
!  to hold R. Thus the amount of memory used, as well as the amount of work
!  involved in computing the preconditioner, depends on mi28_rsize. Setting
!  mi28_rsize > 0 generally leads to a higher quality preconditioner than
!  using mi28_rsize = 0, and choosing mi28_rsize >= mi28_lsize is generally 
!  recommended

!$$     INTEGER :: mi28_rsize = 10

!  which linear least squares solver should we use?
     
     INTEGER :: lls_solver
        
!   overall convergence tolerances. The iteration will terminate when the
!     norm of the gradient of the objective function is smaller than 
!       MAX( %stop_g_absolute, %stop_g_relative * norm of the initial gradient
!     or if the step is less than %stop_s

     REAL ( KIND = wp ) :: stop_g_absolute = tenm5
     REAL ( KIND = wp ) :: stop_g_relative = tenm8
!$$     REAL ( KIND = wp ) :: stop_s = epsmch

!   try to pick a good initial trust-region radius using %advanced_start
!    iterates of a variant on the strategy of Sartenaer SISC 18(6)1990:1788-1803
     
!$$     INTEGER :: advanced_start = 0
     
!   initial value for the trust-region radius (-ve => ||g_0||)
     
     REAL ( KIND = wp ) :: initial_radius = hundred
     
!   maximum permitted trust-region radius

     REAL ( KIND = wp ) :: maximum_radius = ten ** 8

!   a potential iterate will only be accepted if the actual decrease
!    f - f(x_new) is larger than %eta_successful times that predicted
!    by a quadratic model of the decrease. The trust-region radius will be
!    increased if this relative decrease is greater than %eta_very_successful
!    but smaller than %eta_too_successful

     REAL ( KIND = wp ) :: eta_successful = ten ** ( - 8 )
     REAL ( KIND = wp ) :: eta_very_successful = point9
     REAL ( KIND = wp ) :: eta_too_successful = two

!   on very successful iterations, the trust-region radius will be increased by
!    the factor %radius_increase, while if the iteration is unsucceful, the 
!    radius will be decreased by a factor %radius_reduce but no more than
!    %radius_reduce_max

     REAL ( KIND = wp ) :: radius_increase = two
     REAL ( KIND = wp ) :: radius_reduce = half
     REAL ( KIND = wp ) :: radius_reduce_max = sixteenth
       
!   the smallest value the onjective function may take before the problem
!    is marked as unbounded

!$$     REAL ( KIND = wp ) :: obj_unbounded = - epsmch ** ( - 2 )

!   the maximum CPU time allowed (-ve means infinite)
     
!$$     REAL ( KIND = wp ) :: cpu_time_limit = - one

!   the maximum elapsed clock time allowed (-ve means infinite)

!$$     REAL ( KIND = wp ) :: clock_time_limit = - one
       
!   is the Hessian matrix of second derivatives available or is access only
!    via matrix-vector products?

!$$     LOGICAL :: hessian_available = .TRUE.

!   use a direct (factorization) or (preconditioned) iterative method to 
!    find the search direction

!$$     LOGICAL :: subproblem_direct = .FALSE.

!   is a retrospective strategy to be used to update the trust-region radius?

!$$     LOGICAL :: retrospective_trust_region = .FALSE.

!   should the radius be renormalized to account for a change in preconditioner?

!$$     LOGICAL :: renormalize_radius = .FALSE.

!   if %space_critical true, every effort will be made to use as little
!    space as possible. This may result in longer computation time
     
!$$     LOGICAL :: space_critical = .FALSE.
       
!   if %deallocate_error_fatal is true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue

!$$     LOGICAL :: deallocate_error_fatal = .FALSE.

!  all output lines will be prefixed by %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'

!$$     CHARACTER ( LEN = 30 ) :: prefix = '""                            '
     
  END TYPE NLLS_control_type

!  - - - - - - - - - - - - - - - - - - - - - - - 
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - 

  TYPE, PUBLIC :: NLLS_inform_type
     
!  return status. See NLLS_solve for details
     
     INTEGER :: status = 0
     
!  the status of the last attempted allocation/deallocation

!$$     INTEGER :: alloc_status = 0

!  the name of the array for which an allocation/deallocation error ocurred

!$$     CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  the total number of iterations performed
     
!$$     INTEGER :: iter = 0
       
!  the total number of CG iterations performed

!$$     INTEGER :: cg_iter = 0

!  the total number of evaluations of the objection function

!$$     INTEGER :: f_eval = 0

!  the total number of evaluations of the gradient of the objection function

!$$     INTEGER :: g_eval = 0

!  the total number of evaluations of the Hessian of the objection function
     
!$$     INTEGER :: h_eval = 0

!  the maximum number of factorizations in a sub-problem solve

!$$     INTEGER :: factorization_max = 0

!  the return status from the factorization

!$$     INTEGER :: factorization_status = 0

!   the maximum number of entries in the factors

!$$     INTEGER ( KIND = long ) :: max_entries_factors = 0

!  the total integer workspace required for the factorization

!$$     INTEGER :: factorization_integer = - 1

!  the total real workspace required for the factorization

!$$     INTEGER :: factorization_real = - 1

!  the average number of factorizations per sub-problem solve

!$$     REAL ( KIND = wp ) :: factorization_average = zero

!  the value of the objective function at the best estimate of the solution 
!   determined by NLLS_solve

!$$     REAL ( KIND = wp ) :: obj = HUGE( one )

!  the norm of the gradient of the objective function at the best estimate 
!   of the solution determined by NLLS_solve

!$$     REAL ( KIND = wp ) :: norm_g = HUGE( one )

!  the total CPU time spent in the package

!$$     REAL :: cpu_total = 0.0
       
!  the CPU time spent preprocessing the problem

!$$     REAL :: cpu_preprocess = 0.0

!  the CPU time spent analysing the required matrices prior to factorization

!$$     REAL :: cpu_analyse = 0.0

!  the CPU time spent factorizing the required matrices
     
!$$     REAL :: cpu_factorize = 0.0
       
!  the CPU time spent computing the search direction

!$$     REAL :: cpu_solve = 0.0

!  the total clock time spent in the package

!$$     REAL ( KIND = wp ) :: clock_total = 0.0
       
!  the clock time spent preprocessing the problem

!$$     REAL ( KIND = wp ) :: clock_preprocess = 0.0
       
!  the clock time spent analysing the required matrices prior to factorization

!$$     REAL ( KIND = wp ) :: clock_analyse = 0.0
       
!  the clock time spent factorizing the required matrices

!$$     REAL ( KIND = wp ) :: clock_factorize = 0.0
     
!  the clock time spent computing the search direction

!$$     REAL ( KIND = wp ) :: clock_solve = 0.0

  END TYPE NLLS_inform_type

  type params_base_type
     ! deliberately empty
  end type params_base_type
  
  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(:), intent(in)  :: x
       double precision, dimension(:), intent(out) :: f
       class(params_base_type), intent(in) :: params
     end subroutine eval_f_type
  end interface

  abstract interface
     subroutine eval_j_type(status, n, m, x, J, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: J
       class(params_base_type), intent(in) :: params
     end subroutine eval_j_type
  end interface

  
contains


  SUBROUTINE RAL_NLLS( n, m, X,                   & 
                       eval_F, eval_J, params,    &
                       status, options )
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

    USE ISO_FORTRAN_ENV
    INTEGER( int32 ), INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( NLLS_inform_type ), INTENT( OUT ) :: status
    TYPE( NLLS_control_type ), INTENT( IN ) :: options
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    class( params_base_type ) :: params
      
    integer :: jstatus=0, fstatus=0
    integer :: i
    real(wp), DIMENSION(m*n) :: J, Jnew
    real(wp), DIMENSION(m)   :: f, fnew
    real(wp), DIMENSION(n)   :: d, g, Xnew
    real(wp) :: Delta, rho, normJF0, normF0, md

    if ( options%print_level >= 3 )  write( options%out , 3000 ) 

    Delta = options%initial_radius
    
    call eval_J(jstatus, n, m, X, J, params)
    if (jstatus > 0) write( options%out, 1010) jstatus
    call eval_F(fstatus, n, m, X, f, params)
    if (fstatus > 0) write( options%out, 1020) fstatus
    
    normF0 = norm2(f)

!    g = -J^Tf
!    call dgemv('T',m,n,-1.0,J,m,f,1,0.0,g,1)
!    g = - matmul(transpose(J),f)
    call mult_Jt(J,n,m,f,g)
    g = -g
    normJF0 = norm2(g)
    
    main_loop: do i = 1,options%maxit
       
       if ( options%print_level >= 3 )  write( options%out , 3030 ) i
       
       !+++++++++++++++++++++++!
       ! Calculate the step... !
       !+++++++++++++++++++++++!

       call calculate_step(J,f,g,n,m,Delta,d,options)
       
       !++++++++++++++++++!
       ! Accept the step? !
       !++++++++++++++++++!

       Xnew = X + d;
       call eval_J(jstatus, n, m, Xnew, Jnew, params)
       if (jstatus > 0) write( options%out, 1010) jstatus
       call eval_F(fstatus, n, m, Xnew, fnew, params)
       if (fstatus > 0) write( options%out, 1020) fstatus
       
       call evaluate_model(f,J,d,md,m,n,options)

       rho = ( norm2(f)**2 - norm2(fnew)**2 ) / &
             ( norm2(f)**2 - md)
       
       if (rho > options%eta_successful) then
          X = Xnew;
          J = Jnew;
          f = fnew;
!          g = - matmul(transpose(J),f);
          ! g = -J^Tf
!          call dgemv('T',m,n,-1.0,J,m,f,1,0.0,g,1)          
          call mult_Jt(J,n,m,f,g)
          g = -g

          if (options%print_level >=3) write(options%out,3010) 0.5 * norm2(f)**2
          if (options%print_level >=3) write(options%out,3060) norm2(g)/norm2(f)

          !++++++++++++++++++!
          ! Test convergence !
          !++++++++++++++++++!
                   
          if ( norm2(f) <= options%stop_g_absolute + &
               options%stop_g_relative * normF0) then
             if (options%print_level > 0 ) write(options%out,3020) i
             return
          end if
          if ( (norm2(g)/norm2(f)) <= options%stop_g_absolute + &
               options%stop_g_relative * (normJF0/normF0)) then
             if (options%print_level > 0 ) write(options%out,3020) i
             return
          end if
          
       end if
          
       !++++++++++++++++++++++!
       ! Update the TR radius !
       !++++++++++++++++++++++!

       if ( (rho > options%eta_very_successful) &
            .and. (rho < options%eta_too_successful) ) then
          if (options%print_level >=3) write(options%out,3040)
          Delta = max(options%maximum_radius, options%radius_increase * Delta )
       else if (rho < options%eta_successful) then
          if (options%print_level >=3) write(options%out,3050)
          Delta = max( options%radius_reduce, options%radius_reduce_max) * Delta
       end if
       
     end do main_loop

    if (options%print_level > 0 ) write(options%out,1030) 

    RETURN

! Non-executable statements

! print level > 0

1010 FORMAT('Error code from eval_J, status = ',I6)
1020 FORMAT('Error code from eval_J, status = ',I6)
1030 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

! print level > 1

! print level > 2
3000 FORMAT(/,'* Running RAL_NLLS *')
3010 FORMAT('0.5 ||f||^2 = ',ES12.4)
3020 FORMAT('RAL_NLLS converged at iteration ',I6)
3030 FORMAT('== Starting iteration ',i0,' ==')
3040 FORMAT('Very successful step -- increasing Delta')
3050 FORMAT('Successful step -- decreasing Delta')
3060 FORMAT('||J''f||/||f|| = ',ES12.4)
!  End of subroutine RAL_NLLS

  END SUBROUTINE RAL_NLLS
  
  SUBROUTINE calculate_step(J,f,g,n,m,Delta,d,options)
       
! -------------------------------------------------------
! calculate_step, find the next step in the optimization
! -------------------------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: options

     select case (options%nlls_method)
        
     case (1) ! Powell's dogleg
        call dogleg(J,f,g,n,m,Delta,d,options)
        
     case default
        
        if ( options%print_level > 0 ) then
           write(options%error,'(a)') 'Error: unknown value of options%nlls_method'
        end if

     end select

   END SUBROUTINE calculate_step


   SUBROUTINE dogleg(J,f,g,n,m,Delta,d,options)
! -----------------------------------------
! dogleg, implement Powell's dogleg method
! -----------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: options

     real(wp) :: alpha, beta
     real(wp) :: d_sd(n), d_gn(n), ghat(n)
     ! todo: would it be cheaper to allocate this memory in the top loop?
     integer :: slls_status, fb_status
     real(wp) :: Jg(m)
     
!     Jg = 0.0
!     call dgemv('N',m,n,1.0,J,m,g,1,0.0,Jg,1)
     call mult_J(J,n,m,g,Jg)

     alpha = norm2(g)**2 / norm2( Jg )**2
       
     d_sd = alpha * g;

     ! Solve the linear problem...
     select case (options%model)
     case (1)
        ! linear model...
        call solve_LLS(J,f,n,m,options%lls_solver,d_gn,slls_status)
     case default
        if (options%print_level> 0) then
           write(options%error,'(a)') 'Error: model not supported in dogleg'
        end if
        return
     end select

     
     if (norm2(d_gn) <= Delta) then
        d = d_gn
     else if (norm2( alpha * d_sd ) >= Delta) then
        d = (Delta / norm2(d_sd) ) * d_sd
     else
        ghat = d_gn - alpha * d_sd
        call findbeta(d_sd,ghat,alpha,Delta,beta,fb_status)
        d = alpha * d_sd + beta * ghat
     end if
     
   END SUBROUTINE dogleg
     

     SUBROUTINE solve_LLS(J,f,n,m,method,d_gn,status)
       
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------

       REAL(wp), DIMENSION(:), INTENT(IN) :: J
       REAL(wp), DIMENSION(:), INTENT(IN) :: f
       INTEGER, INTENT(IN) :: method, n, m
       REAL(wp), DIMENSION(:), INTENT(OUT) :: d_gn
       INTEGER, INTENT(OUT) :: status

       character(1) :: trans = 'N'
       integer :: nrhs = 1, lwork, lda, ldb
       real(wp), allocatable :: temp(:), work(:)
       REAL(wp), DIMENSION(:,:), allocatable :: Jlls
       integer :: i, jit

       lda = m
       ldb = max(m,n)
       allocate(temp(max(m,n)))
       temp(1:m) = f(1:m)
       lwork = max(1, min(m,n) + max(min(m,n), nrhs)*4)
       allocate(work(lwork))
       allocate(Jlls(m,n))

       do i = 1,n
          do jit = 1,m
             Jlls(jit,i) = J((i-1) * m + jit) ! We need to take a copy as dgels overwrites J
          end do
       end do

       call dgels(trans, m, n, nrhs, Jlls, lda, temp, ldb, work, lwork, status)

       d_gn = -temp(1:n)
              
     END SUBROUTINE solve_LLS
     
     SUBROUTINE findbeta(d_sd,ghat,alpha,Delta,beta,status)

!  -----------------------------------------------------------------
!  findbeta, a subroutine to find the optimal beta such that 
!   || d || = Delta
!  -----------------------------------------------------------------

     real(wp), dimension(:), intent(in) :: d_sd, ghat
     real(wp), intent(in) :: alpha, Delta
     real(wp), intent(out) :: beta
     integer, intent(out) :: status
     
     real(wp) :: a, b, c, discriminant

     a = norm2(ghat)**2
     b = 2.0 * alpha * dot_product( ghat, d_sd)
     c = ( alpha * norm2( d_sd ) )**2 - Delta
     
     discriminant = b**2 - 4 * a * c
     if ( discriminant < 0) then
        status = 1
        return
     else
        beta = (-b + sqrt(discriminant)) / (2.0 * a)
     end if

     END SUBROUTINE findbeta

     
     subroutine evaluate_model(f,J,d,md,m,n,options)
! --------------------------------------------------
! evaluate_model, evaluate the model at the point d
! --------------------------------------------------       

       real(wp), intent(in)  :: f(:), d(:), J(:)
       integer :: m,n
       real(wp), intent(out) :: md
       TYPE( NLLS_control_type ), INTENT( IN ) :: options

       real(wp) :: Jd(m)
       
       !Jd = J*d
!       call dgemv('N',m,n,1.0,J,m,d,1,0.0,Jd,1)
       call mult_J(J,n,m,d,Jd)

       select case (options%model)
       case (1) ! first-order (no Hessian)
          md = norm2(f + Jd)**2
       case (2) ! second-order (exact hessian)
          if (options%print_level > 0) then
             write(options%error,'(a)') 'Second order methods to be implemented'
             return
          end if
       case (3) ! barely second-order (identity Hessian)
          md = norm2(f + Jd + dot_product(d,d))**2
       case default
          if (options%print_level > 0) then
             write(options%error,'(a)') 'Error: unsupported model specified'
          end if
          return
          
       end select

     end subroutine evaluate_model


     subroutine mult_J(J,n,m,x,Jx)
       real(wp), intent(in) :: J(:), x(:)
       integer, intent(in) :: n,m
       real(wp), intent(out) :: Jx(:)

       real(wp), allocatable :: Jarray(:,:)
       integer :: ii, jj
       
       allocate(Jarray(m,n))
       do ii = 1,n
          do jj = 1,m
             Jarray(jj,ii) = J((ii-1) * m  + jj)
          end do
       end do

       Jx = matmul(Jarray,x)
       
     end subroutine mult_J

     subroutine mult_Jt(J,n,m,x,Jtx)
       real(wp), intent(in) :: J(:), x(:)
       integer, intent(in) :: n,m
       real(wp), intent(out) :: Jtx(:)

       real(wp), allocatable :: Jtarray(:,:)
       integer :: ii, jj
       
       allocate(Jtarray(n,m))

       do ii = 1,n
          do jj = 1,m
             Jtarray(ii,jj) = J((ii-1) * m  + jj)
          end do
       end do

       Jtx = matmul(Jtarray,x)
       
     end subroutine mult_Jt

end module nlls_module
