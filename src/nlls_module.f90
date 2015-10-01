! nlls_module :: a nonlinear least squares solver

module nlls_module

  implicit none

  integer, parameter :: wp = kind(1.0d0)
  
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

     INTEGER :: start_print = - 1

!   any printing will stop on this iteration

     INTEGER :: stop_print = - 1

!   the number of iterations between printing

     INTEGER :: print_gap = 1

!   the maximum number of iterations performed

     INTEGER :: maxit = 100

!   removal of the file alive_file from unit alive_unit terminates execution

     INTEGER :: alive_unit = 40
     CHARACTER ( LEN = 30 ) :: alive_file = 'ALIVE.d'

!   non-monotone <= 0 monotone strategy used, anything else non-monotone
!     strategy with this history length used

     INTEGER :: non_monotone = 1

!   specify the model used. Possible values are
!
!      0  dynamic (*not yet implemented*)
!      1  first-order (no Hessian)
!      2  second-order (exact Hessian)
!      3  barely second-order (identity Hessian)
!      4  secant second-order (sparsity-based)
!      5  secant second-order (limited-memory BFGS, with %lbfgs_vectors history)
!      6  secant second-order (limited-memory SR1, with %lbfgs_vectors history)

     INTEGER :: model = 2

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

     INTEGER :: norm = 1

!   specify the semi-bandwidth of the band matrix P if required

     INTEGER :: semi_bandwidth = 5

!   number of vectors used by the L-BFGS matrix P if required

     INTEGER :: lbfgs_vectors = 10

!   number of vectors used by the sparsity-based secant Hessian if required

     INTEGER :: max_dxg = 100

!   number of vectors used by the Lin-More' incomplete factorization 
!    matrix P if required

     INTEGER :: icfs_vectors = 10

!  the maximum number of fill entries within each column of the incomplete 
!  factor L computed by HSL_MI28. In general, increasing mi28_lsize improves
!  the quality of the preconditioner but increases the time to compute
!  and then apply the preconditioner. Values less than 0 are treated as 0

     INTEGER :: mi28_lsize = 10

!  the maximum number of entries within each column of the strictly lower 
!  triangular matrix R used in the computation of the preconditioner by 
!  HSL_MI28.  Rank-1 arrays of size mi28_rsize *  n are allocated internally 
!  to hold R. Thus the amount of memory used, as well as the amount of work
!  involved in computing the preconditioner, depends on mi28_rsize. Setting
!  mi28_rsize > 0 generally leads to a higher quality preconditioner than
!  using mi28_rsize = 0, and choosing mi28_rsize >= mi28_lsize is generally 
!  recommended

     INTEGER :: mi28_rsize = 10
        
!   overall convergence tolerances. The iteration will terminate when the
!     norm of the gradient of the objective function is smaller than 
!       MAX( %stop_g_absolute, %stop_g_relative * norm of the initial gradient
!     or if the step is less than %stop_s

     REAL ( KIND = wp ) :: stop_g_absolute = tenm5
     REAL ( KIND = wp ) :: stop_g_relative = tenm8
     REAL ( KIND = wp ) :: stop_s = epsmch

!   try to pick a good initial trust-region radius using %advanced_start
!    iterates of a variant on the strategy of Sartenaer SISC 18(6)1990:1788-1803
     
     INTEGER :: advanced_start = 0
     
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

     REAL ( KIND = wp ) :: obj_unbounded = - epsmch ** ( - 2 )

!   the maximum CPU time allowed (-ve means infinite)
     
     REAL ( KIND = wp ) :: cpu_time_limit = - one

!   the maximum elapsed clock time allowed (-ve means infinite)

     REAL ( KIND = wp ) :: clock_time_limit = - one
       
!   is the Hessian matrix of second derivatives available or is access only
!    via matrix-vector products?

     LOGICAL :: hessian_available = .TRUE.

!   use a direct (factorization) or (preconditioned) iterative method to 
!    find the search direction

     LOGICAL :: subproblem_direct = .FALSE.

!   is a retrospective strategy to be used to update the trust-region radius?

     LOGICAL :: retrospective_trust_region = .FALSE.

!   should the radius be renormalized to account for a change in preconditioner?

     LOGICAL :: renormalize_radius = .FALSE.

!   if %space_critical true, every effort will be made to use as little
!    space as possible. This may result in longer computation time
     
     LOGICAL :: space_critical = .FALSE.
       
!   if %deallocate_error_fatal is true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue

     LOGICAL :: deallocate_error_fatal = .FALSE.

!  all output lines will be prefixed by %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'

     CHARACTER ( LEN = 30 ) :: prefix = '""                            '
     
  END TYPE NLLS_control_type

!  - - - - - - - - - - - - - - - - - - - - - - - 
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - 

  TYPE, PUBLIC :: NLLS_inform_type
     
!  return status. See NLLS_solve for details
     
     INTEGER :: status = 0
     
!  the status of the last attempted allocation/deallocation

     INTEGER :: alloc_status = 0

!  the name of the array for which an allocation/deallocation error ocurred

     CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  the total number of iterations performed
     
     INTEGER :: iter = 0
       
!  the total number of CG iterations performed

     INTEGER :: cg_iter = 0

!  the total number of evaluations of the objection function

     INTEGER :: f_eval = 0

!  the total number of evaluations of the gradient of the objection function

     INTEGER :: g_eval = 0

!  the total number of evaluations of the Hessian of the objection function
     
     INTEGER :: h_eval = 0

!  the maximum number of factorizations in a sub-problem solve

     INTEGER :: factorization_max = 0

!  the return status from the factorization

     INTEGER :: factorization_status = 0

!   the maximum number of entries in the factors

     INTEGER ( KIND = long ) :: max_entries_factors = 0

!  the total integer workspace required for the factorization

     INTEGER :: factorization_integer = - 1

!  the total real workspace required for the factorization

     INTEGER :: factorization_real = - 1

!  the average number of factorizations per sub-problem solve

     REAL ( KIND = wp ) :: factorization_average = zero

!  the value of the objective function at the best estimate of the solution 
!   determined by NLLS_solve

     REAL ( KIND = wp ) :: obj = HUGE( one )

!  the norm of the gradient of the objective function at the best estimate 
!   of the solution determined by NLLS_solve

     REAL ( KIND = wp ) :: norm_g = HUGE( one )

!  the total CPU time spent in the package

     REAL :: cpu_total = 0.0
       
!  the CPU time spent preprocessing the problem

     REAL :: cpu_preprocess = 0.0

!  the CPU time spent analysing the required matrices prior to factorization

     REAL :: cpu_analyse = 0.0

!  the CPU time spent factorizing the required matrices
     
     REAL :: cpu_factorize = 0.0
       
!  the CPU time spent computing the search direction

     REAL :: cpu_solve = 0.0

!  the total clock time spent in the package

     REAL ( KIND = wp ) :: clock_total = 0.0
       
!  the clock time spent preprocessing the problem

     REAL ( KIND = wp ) :: clock_preprocess = 0.0
       
!  the clock time spent analysing the required matrices prior to factorization

     REAL ( KIND = wp ) :: clock_analyse = 0.0
       
!  the clock time spent factorizing the required matrices

     REAL ( KIND = wp ) :: clock_factorize = 0.0
     
!  the clock time spent computing the search direction

     REAL ( KIND = wp ) :: clock_solve = 0.0

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
    REAL( real64 ), DIMENSION( n ), INTENT( INOUT ) :: X
    INTEGER( int32 ), INTENT( OUT ) :: status
     
!  Interface blocks (e.g.)

    INTERFACE
!      SUBROUTINE eval_F( status, X, userdata, f )
       SUBROUTINE eval_F( status, X, f )
         USE ISO_FORTRAN_ENV
         !      USE GALAHAD_NLPT_double, ONLY: NLPT_userdata_type
         INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
         INTEGER ( int32 ), INTENT( OUT ) :: status
         REAL ( real64 ), INTENT( OUT ) :: f
         REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X
         !      TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
       END SUBROUTINE eval_F
    END INTERFACE

     INTERFACE
!      SUBROUTINE eval_J( status, X, userdata, J )
        SUBROUTINE eval_J( status, X, J )
          USE ISO_FORTRAN_ENV
          !      USE GALAHAD_NLPT_double, ONLY: NLPT_userdata_type
          INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
          INTEGER ( int32 ), INTENT( OUT ) :: status
          REAL ( real64 ), INTENT( OUT ) :: f
          REAL ( real64 ), DIMENSION( : , : ),INTENT( IN ) :: J
          !      TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
        END SUBROUTINE eval_J
     END INTERFACE

     RETURN

!  End of subroutine RAL_NLLS

     END SUBROUTINE RAL_NLLS


end module nlls_module
