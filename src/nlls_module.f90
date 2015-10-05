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

!$$     INTEGER :: model = 2

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

!   specify the method used to solve the nlls problem
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

contains


  SUBROUTINE RAL_NLLS( n, m, X, Work_int, len_work_int,                     &
                       Work_real, len_work_real,                            &
                       eval_F, eval_J,                                      &
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
    INTEGER( int32 ), INTENT( IN ) :: n, m, len_work_int, len_work_real
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    INTEGER( int32), INTENT( IN ) :: Work_int(len_work_int)
    REAL( wp ), INTENT( IN ) :: Work_real(len_work_real)
    TYPE( NLLS_inform_type ), INTENT( OUT ) :: status
    TYPE( NLLS_control_type ), INTENT( IN ) :: options

!  Interface blocks (e.g.)

    INTERFACE
       SUBROUTINE eval_F( status, X, f )
         USE ISO_FORTRAN_ENV
         
         INTEGER ( int32 ), INTENT( OUT ) :: status
         REAL ( real64 ), DIMENSION( : ),INTENT( OUT ) :: f
         REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X
         
       END SUBROUTINE eval_F
    END INTERFACE

    INTERFACE
       SUBROUTINE eval_J( status, X, J )
         USE ISO_FORTRAN_ENV

         INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
         INTEGER ( int32 ), INTENT( OUT ) :: status
         REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X
         REAL ( real64 ), DIMENSION( : , : ),INTENT( OUT ) :: J
       END SUBROUTINE eval_J
    END INTERFACE
    
    integer :: jstatus=0, fstatus=0
    integer :: i
    real(wp), DIMENSION(m,n) :: J, Jnew
    real(wp), DIMENSION(m) :: f, fnew
    real(wp), DIMENSION(n) :: d, g, Xnew
    real(wp) :: Delta, rho, normJF0, normF0

    if ( options%print_level >= 3 )  write( options%out , 3000 ) 

    Delta = options%initial_radius
    
    call eval_J(jstatus, X, J)
    if (jstatus > 0) write( options%out, 1010) jstatus
    call eval_F(fstatus, X, f)
    if (fstatus > 0) write( options%out, 1020) fstatus
    
    normF0 = norm2(f)
    g = - matmul(transpose(J),f);
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
       call eval_J(jstatus, Xnew, Jnew)
       if (jstatus > 0) write( options%out, 1010) jstatus
       call eval_F(fstatus, Xnew, fnew)
       if (fstatus > 0) write( options%out, 1020) fstatus
       
       rho = ( norm2(f)**2 - norm2(fnew)**2 ) / &
             ( norm2(f)**2 - norm2(f + matmul(J,d))**2)
       
       if (rho > options%eta_successful) then
          X = Xnew;
          J = Jnew;
          f = fnew;
          g = - matmul(transpose(J),f);
          
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

     REAL(wp), intent(in) :: J(:,:), f(:), g(:), Delta
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

     REAL(wp), intent(in) :: J(:,:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: options

     real(wp) :: alpha, beta
     real(wp) :: d_sd(n), d_gn(n), ghat(n)
     ! todo: would it be cheaper to allocate this memory in the top loop?
     integer :: slls_status, fb_status

     alpha = norm2(g)**2 / norm2( matmul(J,g) )**2
       
     d_sd = alpha * g;
     call solve_LLS(J,f,n,m,options%lls_solver,d_gn,slls_status)
     
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

       REAL(wp), DIMENSION(:,:), INTENT(IN) :: J
       REAL(wp), DIMENSION(:), INTENT(IN) :: f
       INTEGER, INTENT(IN) :: method, n, m
       REAL(wp), DIMENSION(:), INTENT(OUT) :: d_gn
       INTEGER, INTENT(OUT) :: status

       character(1) :: trans = 'N'
       integer :: nrhs = 1, lwork, lda, ldb
       real(wp), allocatable :: temp(:), work(:)
       REAL(wp), DIMENSION(:,:), allocatable :: Jlls
       integer :: i

       lda = m
       ldb = max(m,n)
       allocate(temp(max(m,n)))
       temp(1:m) = f(1:m)
       lwork = max(1, min(m,n) + max(min(m,n), nrhs)*4)
       allocate(work(lwork))
       
       Jlls = J ! We need to take a copy as dgels overwrites J
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

     
     SUBROUTINE eval_F( status, X, f )

!  -------------------------------------------------------------------
!  eval_F, a subroutine for evaluating the function f at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER ( int32 ), INTENT( OUT ) :: status
       REAL ( real64 ), DIMENSION( : ),INTENT( OUT ) :: f
       REAL ( real64 ), DIMENSION( : ),INTENT( IN )  :: X

! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer, parameter :: num_observations = 67
       integer :: i
       real( wp ) :: x_data( num_observations )
       real( wp ) :: y_data( num_observations )

       call generate_data_example(x_data,y_data,num_observations)

! then, let's work this into the format we need
! X(1) = m, X(2) = c
       do i = 1,num_observations
          f(i) = y_data(i) - exp( X(1) * x_data(i) + X(2) )
       end do
       

!!$! let's use Powell's function for now....
!!$       f(1) = X(1) + 10.0 * X(2)
!!$       f(2) = sqrt(5.0) * (X(3) - X(4))
!!$       f(3) = ( X(2) - 2.0 * X(3) )**2
!!$       f(4) = sqrt(10.0) * ( X(1) - X(4) )**2
       
! end of subroutine eval_F
       
     END SUBROUTINE eval_F

     SUBROUTINE eval_J( status, X, J )

!  -------------------------------------------------------------------
!  eval_J, a subroutine for evaluating the Jacobian J at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER ( int32 ), INTENT( OUT ) :: status
       REAL ( real64 ), DIMENSION( : , : ),INTENT( OUT ) :: J
       REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X

! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer, parameter :: num_observations = 67
       integer :: i
       real( wp ) :: x_data( num_observations ), y_data( num_observations )

! First, let's re-generate the data....
       call generate_data_example(x_data,y_data,num_observations)

! then, let's work this into the format we need
! X(1) = m, X(2) = c
       do i = 1,num_observations
          J(i,1) =  - x_data(i) * exp( X(1) * x_data(i) + X(2) )
          J(i,2) =  - exp( X(1) * x_data(i) + X(2) )
       end do
       

! end of subroutine eval_J

!!$       ! initialize to zeros...
!!$       J(1:4,1:4) = 0.0
!!$       
!!$       ! enter non-zeros values
!!$       J(1,1) = 1.0
!!$       J(1,2) = 10.0
!!$       J(2,3) = sqrt(5.0)
!!$       J(2,4) = -sqrt(5.0)
!!$       J(3,2) = 2.0 * (X(2) - 2.0 * X(3))
!!$       J(3,3) = -4.0 * (X(2) - 2.0 * X(3)) 
!!$       J(4,1) = sqrt(10.0) * 2.0 * (X(1) - X(4))
!!$       J(4,4) = - sqrt(10.0) * 2.0 * (X(1) - X(4))

     END SUBROUTINE eval_J

     subroutine generate_data_example(x_data,y_data,num_observations)

       real(wp), intent(out) :: x_data(:), y_data(:)
       integer, intent(in) :: num_observations

       ! First, let's get the data
       ! Generated with the code
       !  randn('seed', 23497);
       !   m = 0.3;
       !   c = 0.1;
       !   x_data = [0:0.075:5];
       !   y = exp(m * x_data + c);
       !   noise = randn(size(x_data)) * 0.2;
       !   y_data = y + noise;
       ! (c.f. https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/curve_fitting.cc)
       x_data = (/ 0.0, &
   0.075000000000000, &
   0.150000000000000, &
   0.225000000000000, &
   0.300000000000000, &
   0.375000000000000, &
   0.450000000000000, &
   0.525000000000000, &
   0.600000000000000, &
   0.675000000000000, &
   0.750000000000000, &
   0.825000000000000, &
   0.900000000000000, &
   0.975000000000000, &
   1.050000000000000, &
   1.125000000000000, &
   1.200000000000000, &
   1.275000000000000, &
   1.350000000000000, &
   1.425000000000000, &
   1.500000000000000, &
   1.575000000000000, &
   1.650000000000000, &
   1.725000000000000, &
   1.800000000000000, &
   1.875000000000000, &
   1.950000000000000, &
   2.025000000000000, &
   2.100000000000000, &
   2.175000000000000, &
   2.250000000000000, &
   2.325000000000000, &
   2.400000000000000, &
   2.475000000000000, &
   2.550000000000000, &
   2.625000000000000, &
   2.700000000000000, &
   2.775000000000000, &
   2.850000000000000, &
   2.925000000000000, &
   3.000000000000000, &
   3.075000000000000, &
   3.150000000000000, &
   3.225000000000001, &
   3.300000000000000, &
   3.375000000000000, &
   3.450000000000000, &
   3.525000000000000, &
   3.600000000000001, &
   3.675000000000000, &
   3.750000000000000, &
   3.825000000000000, &
   3.900000000000000, &
   3.975000000000000, &
   4.050000000000001, &
   4.125000000000000, &
   4.200000000000000, &
   4.275000000000000, &
   4.350000000000001, &
   4.425000000000000, &
   4.500000000000000, &
   4.575000000000000, &
   4.650000000000000, &
   4.725000000000001, &
   4.800000000000000, &
   4.875000000000000, &
   4.950000000000000 /)

       y_data = (/ 0.907946872110432, &
   1.199579396036134, &
   1.060092431384317, &
   1.298370500472354, &
   0.952768858414788, &
   1.209665290655204, &
   1.256912538155493, &
   1.163922146095987, &
   1.004877938808100, &
   1.205944250961060, &
   0.952693297695969, &
   1.449662692280761, &
   1.402015259144406, &
   1.378094012325746, &
   1.560882147577552, &
   1.437185539058121, &
   1.559853079888265, &
   1.877814947316832, &
   1.818781749024682, &
   1.375546045112591, &
   1.233967904388409, &
   1.887793124397751, &
   1.610237096463521, &
   1.787032484792262, &
   1.850015127982676, &
   2.120553361509177, &
   1.942913663511919, &
   2.106517132599766, &
   2.271787117356578, &
   1.727554346001754, &
   2.002909500898113, &
   1.975837413903495, &
   2.337446525801909, &
   1.960190841677278, &
   2.447097025572309, &
   2.161663720225506, &
   2.748798529374621, &
   2.507814238594416, &
   2.423769408403069, &
   2.578119353028746, &
   2.460310096221557, &
   2.638362783992324, &
   2.765540456237868, &
   2.837165966564409, &
   3.179711963042789, &
   3.245315453091675, &
   3.289631922410174, &
   3.360995198615834, &
   3.470489725998371, &
   3.169513520153466, &
   3.363740517933189, &
   3.665288099084969, &
   3.620334359722351, &
   4.018911445550667, &
   3.512715166706162, &
   3.874661411575566, &
   4.197746303653517, &
   3.703511523106007, &
   4.076351488309604, &
   4.056340365649961, &
   4.297751562451419, &
   4.373076571153739, &
   4.577093065941748, &
   4.856619059058190, &
   4.927350280596274, &
   4.703122139742729, &
   4.870205182453842 /)
       
     end subroutine generate_data_example


end module nlls_module
