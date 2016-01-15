! nlls_module :: a nonlinear least squares solver

module nlls_module

  use RAL_NLLS_DTRS_double

  implicit none

  integer, parameter :: wp = kind(1.0d0)
  integer, parameter :: long = selected_int_kind(8)
  real (kind = wp), parameter :: tenm3 = 1.0e-3
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! P R I N T I N G   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! M A I N   R O U T I N E   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!      1  Gauss-Newton (no 2nd derivatives)
!      2  second-order (exact Hessian)
!      3  barely second-order (identity Hessian)
!      4  secant second-order (sparsity-based)
!      5  secant second-order (limited-memory BFGS, with %lbfgs_vectors history)
!      6  secant second-order (limited-memory SR1, with %lbfgs_vectors history)
!      7  hybrid (Gauss-Newton until gradient small, then Newton)
!      8  hybrid (using Madsen, Nielsen and Tingleff's method)    
 
     INTEGER :: model = 9

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
!
!$$     INTEGER :: norm = 1


     INTEGER :: nlls_method = 1

!   specify the method used to solve the trust-region sub problem
!      1 Powell's dogleg
!      2 AINT method (of Yuji Nat.)
!      3 More-Sorensen
!      4 Galahad's DTRS
!      ...

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
     
!   should we scale the initial trust region radius?
     
     integer :: relative_tr_radius = 0

!   if relative_tr_radius == 1, then pick a scaling parameter
!   Madsen, Nielsen and Tingleff say pick this to be 1e-6, say, if x_0 is good,
!   otherwise 1e-3 or even 1 would be good starts...
     
     real (kind = wp) :: initial_radius_scale = 1.0!tenm3

!   if relative_tr_radius /= 1, then set the 
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
!    the factor %radius_increase, while if the iteration is unsuccessful, the 
!    radius will be decreased by a factor %radius_reduce but no more than
!    %radius_reduce_max

     REAL ( KIND = wp ) :: radius_increase = two
     REAL ( KIND = wp ) :: radius_reduce = half
     REAL ( KIND = wp ) :: radius_reduce_max = sixteenth
       
!   the smallest value the objective function may take before the problem
!    is marked as unbounded

!$$     REAL ( KIND = wp ) :: obj_unbounded = - epsmch ** ( - 2 )

!   if model=7, then the value with which we switch on second derivatives
     
     real ( kind = wp ) :: hybrid_switch = 0.1_wp

!   the maximum CPU time allowed (-ve means infinite)
     
!$$     REAL ( KIND = wp ) :: cpu_time_limit = - one

!   the maximum elapsed clock time allowed (-ve means infinite)

!$$     REAL ( KIND = wp ) :: clock_time_limit = - one
 
!   shall we use explicit second derivatives, or approximate using a secant 
!   method
     
     LOGICAL :: exact_second_derivatives = .true.
      
!   is the Hessian matrix of second derivatives available or is access only
!    via matrix-vector products?

!     LOGICAL :: hessian_available = .TRUE.

!   use a direct (factorization) or (preconditioned) iterative method to 
!    find the search direction

!$$     LOGICAL :: subproblem_direct = .FALSE.

!   use a factorization (dsyev) to find the smallest eigenvalue for the subproblem
!    solve? (alternative is an iterative method (dsyevx)
     LOGICAL :: subproblem_eig_fact = .FALSE. ! undocumented....
     

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! M O R E - S O R E N S E N   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

    integer  :: more_sorensen_maxits = 500
    real(wp) :: more_sorensen_shift = 1e-13
    real(wp) :: more_sorensen_tiny = 10.0 * epsmch
    real(wp) :: more_sorensen_tol = 1e-3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! H Y B R I D   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what's the tolerance such that ||J^T f || < tol * 0.5 ||f||^2 triggers a switch
    real(wp) :: hybrid_tol = 0.02

! how many successive iterations does the above condition need to hold before we switch?
    integer  :: hybrid_switch_its = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! O U T P U T   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Shall we output progess vectors at termination of the routine?
     logical :: output_progress_vectors = .false.

  END TYPE NLLS_control_type

!  - - - - - - - - - - - - - - - - - - - - - - - 
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - 

  TYPE, PUBLIC :: NLLS_inform_type
     
!  return status
!   1 -- maximum number of iterations reached
!   2 -- error from evaluating a function/Jacobian/Hessian
!   3 -- unsupported choice of model
!   4 -- error return from an lapack routine
     
     INTEGER :: status = 0
     
!  the status of the last attempted allocation/deallocation

     INTEGER :: alloc_status = 0

!  the name of the array for which an allocation/deallocation error ocurred

!$$     CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  the total number of iterations performed
     
     INTEGER :: iter
       
!  the total number of CG iterations performed

!$$     INTEGER :: cg_iter = 0

!  the total number of evaluations of the objective function

     INTEGER :: f_eval = 0

!  the total number of evaluations of the gradient of the objective function

     INTEGER :: g_eval = 0

!  the total number of evaluations of the Hessian of the objective function
     
     INTEGER :: h_eval = 0

!  test on the size of f satisfied?
     
     integer :: convergence_normf = 0

!  test on the size of the gradient satisfied?
     
     integer :: convergence_normg = 0
     
!  vector of residuals 
     
     real(wp), allocatable :: resvec(:)

!  vector of gradients 
     
     real(wp), allocatable :: gradvec(:)

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

     REAL ( KIND = wp ) :: obj = HUGE( one )

!  the norm of the gradient of the objective function at the best estimate 
!   of the solution determined by NLLS_solve

     REAL ( KIND = wp ) :: norm_g = HUGE( one )

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
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: f
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

  abstract interface
     subroutine eval_hf_type(status, n, m, x, f, h, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: f
       double precision, dimension(*), intent(out) :: h
       class(params_base_type), intent(in) :: params
     end subroutine eval_hf_type
  end interface


  ! define types for workspace arrays.
    
    type, private :: max_eig_work ! workspace for subroutine max_eig
       real(wp), allocatable :: alphaR(:), alphaI(:), beta(:), vr(:,:)
       real(wp), allocatable :: work(:), ew_array(:)
       integer, allocatable :: nullindex(:)
       logical, allocatable :: vecisreal(:)
    end type max_eig_work

    type, private :: solve_general_work ! workspace for subroutine solve_general
       real(wp), allocatable :: A(:,:)
       integer, allocatable :: ipiv(:)
    end type solve_general_work

    type, private :: evaluate_model_work ! workspace for subroutine evaluate_model
       real(wp), allocatable :: Jd(:), Hd(:)
    end type evaluate_model_work

    type, private :: solve_LLS_work ! workspace for subroutine solve_LLS
       real(wp), allocatable :: temp(:), work(:), Jlls(:)
    end type solve_LLS_work
    
    type, private :: min_eig_symm_work ! workspace for subroutine min_eig_work
       real(wp), allocatable :: A(:,:), work(:), ew(:)
       integer, allocatable :: iwork(:), ifail(:)
    end type min_eig_symm_work
    
    type, private :: all_eig_symm_work ! workspace for subroutine all_eig_symm
       real(wp), allocatable :: work(:)
    end type all_eig_symm_work
        
    type, private :: solve_dtrs_work ! workspace for subroutine dtrs_work
       real(wp), allocatable :: A(:,:), ev(:,:), ew(:), v(:), v_trans(:), d_trans(:)
       type( all_eig_symm_work ) :: all_eig_symm_ws
    end type solve_dtrs_work
        
    type, private :: more_sorensen_work ! workspace for subroutine more_sorensen
 !      type( solve_spd_work ) :: solve_spd_ws
       real(wp), allocatable :: A(:,:), LtL(:,:), AplusSigma(:,:)
       real(wp), allocatable :: v(:), q(:), y1(:)
!       type( solve_general_work ) :: solve_general_ws
       type( min_eig_symm_work ) :: min_eig_symm_ws
    end type more_sorensen_work

    type, private :: AINT_tr_work ! workspace for subroutine AINT_tr
       type( max_eig_work ) :: max_eig_ws
       type( evaluate_model_work ) :: evaluate_model_ws
       type( solve_general_work ) :: solve_general_ws
!       type( solve_spd_work ) :: solve_spd_ws
       REAL(wp), allocatable :: A(:,:), LtL(:,:), v(:), B(:,:), p0(:), p1(:)
       REAL(wp), allocatable :: M0(:,:), M1(:,:), y(:), gtg(:,:), q(:)
    end type AINT_tr_work

    type, private :: dogleg_work ! workspace for subroutine dogleg
       type( solve_LLS_work ) :: solve_LLS_ws
       type( evaluate_model_work ) :: evaluate_model_ws
       real(wp), allocatable :: d_sd(:), d_gn(:), ghat(:), Jg(:)
    end type dogleg_work

    type, private :: calculate_step_work ! workspace for subroutine calculate_step
       type( AINT_tr_work ) :: AINT_tr_ws
       type( dogleg_work ) :: dogleg_ws
       type( more_sorensen_work ) :: more_sorensen_ws
       type( solve_dtrs_work ) :: solve_dtrs_ws
    end type calculate_step_work

    type :: NLLS_workspace ! all workspaces called from the top level
       integer :: first_call = 1
       integer :: iter = 0 
       real(wp) :: normF0, normJF0
       real(wp) :: normJFold, normJF_Newton
       real(wp) :: Delta
       logical :: use_second_derivatives = .false.
       integer :: hybrid_count = 0
       real(wp) :: hybrid_tol = 1.0
       real(wp), allocatable :: fNewton(:), JNewton(:), XNewton(:)
       real(wp), allocatable :: J(:)
       real(wp), allocatable :: f(:), fnew(:)
       real(wp), allocatable :: hf(:)
       real(wp), allocatable :: d(:), g(:), Xnew(:)
       real(wp), allocatable :: y(:), y_sharp(:), g_old(:), g_mixed(:)
       real(wp), allocatable :: resvec(:), gradvec(:)
       type ( calculate_step_work ) :: calculate_step_ws
       type ( evaluate_model_work ) :: evaluate_model_ws
    end type NLLS_workspace

contains


  SUBROUTINE RAL_NLLS( n, m, X,                   & 
                       eval_F, eval_J, eval_HF,   & 
                       params,                    &
                       info, control )
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( NLLS_inform_type ), INTENT( OUT ) :: info
    TYPE( NLLS_control_type ), INTENT( IN ) :: control
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
      
    integer  :: i
    
    type ( NLLS_workspace ) :: w
    
!!$    write(*,*) 'Controls in:'
!!$    write(*,*) control
!!$    write(*,*) 'error = ',control%error
!!$    write(*,*) 'out = ', control%out
!!$    write(*,*) 'print_level = ', control%print_level
!!$    write(*,*) 'maxit = ', control%maxit
!!$    write(*,*) 'model = ', control%model
!!$    write(*,*) 'nlls_method = ', control%nlls_method
!!$    write(*,*) 'lls_solver = ', control%lls_solver
!!$    write(*,*) 'stop_g_absolute = ', control%stop_g_absolute
!!$    write(*,*) 'stop_g_relative = ', control%stop_g_relative     
!!$    write(*,*) 'initial_radius = ', control%initial_radius
!!$    write(*,*) 'maximum_radius = ', control%maximum_radius
!!$    write(*,*) 'eta_successful = ', control%eta_successful
!!$    write(*,*) 'eta_very_successful = ',control%eta_very_successful
!!$    write(*,*) 'eta_too_successful = ',control%eta_too_successful
!!$    write(*,*) 'radius_increase = ',control%radius_increase
!!$    write(*,*) 'radius_reduce = ',control%radius_reduce
!!$    write(*,*) 'radius_reduce_max = ',control%radius_reduce_max
!!$    write(*,*) 'hybrid_switch = ',control%hybrid_switch
!!$    write(*,*) 'subproblem_eig_fact = ',control%subproblem_eig_fact
!!$    write(*,*) 'more_sorensen_maxits = ',control%more_sorensen_maxits
!!$    write(*,*) 'more_sorensen_shift = ',control%more_sorensen_shift
!!$    write(*,*) 'more_sorensen_tiny = ',control%more_sorensen_tiny
!!$    write(*,*) 'more_sorensen_tol = ',control%more_sorensen_tol
!!$    write(*,*) 'hybrid_tol = ', control%hybrid_tol
!!$    write(*,*) 'hybrid_switch_its = ', control%hybrid_switch_its
!!$    write(*,*) 'output_progress_vectors = ',control%output_progress_vectors

    main_loop: do i = 1,control%maxit
       
       call ral_nlls_iterate(n, m, X,                   & 
                             eval_F, eval_J, eval_HF,   & 
                             params,                    &
                             info, control, w )
       ! test the returns to see if we've converged
       if (info%status .ne. 0) then 
          return 
       elseif ((info%convergence_normf == 1).or.(info%convergence_normg == 1)) then
          return
       end if
       
     end do main_loop
    
     ! If we reach here, then we're over maxits     
     if (control%print_level > 0 ) write(control%out,1040) 
     info%status = -1
    
     RETURN

! Non-executable statements

! print level > 0

1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

!  End of subroutine RAL_NLLS

  END SUBROUTINE RAL_NLLS
  
  subroutine ral_nlls_iterate(n, m, X,                   & 
                              eval_F, eval_J, eval_HF,   & 
                              params,                    &
                              info, control, w )

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( NLLS_inform_type ), INTENT( OUT ) :: info
    TYPE( NLLS_control_type ), INTENT( IN ) :: control
    type( NLLS_workspace ), INTENT( INOUT ) :: w
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
      
    integer :: jstatus=0, fstatus=0, hfstatus=0, astatus = 0, svdstatus = 0
    integer :: i, no_reductions, max_tr_decrease = 100
    real(wp) :: rho, normJF, normF, normFnew, md, Jmax, JtJdiag
    real(wp) :: FunctionValue, hybrid_tol
    logical :: success, calculate_svd_J
    real(wp) :: s1, sn
    
    ! todo: make max_tr_decrease a control variable

    ! Perform a single iteration of the RAL_NLLS loop
    
    calculate_svd_J = .true. ! todo :: make a control variable 

    if (w%first_call == 1) then
       ! This is the first call...allocate arrays, and get initial 
       ! function evaluations
       if ( control%print_level >= 3 )  write( control%out , 3000 ) 
       ! first, check if n < m
       if (n > m) goto 4070

       call setup_workspaces(w,n,m,control,info%alloc_status)
       if ( info%alloc_status > 0) goto 4000

       call eval_F(fstatus, n, m, X, w%f, params)
       info%f_eval = info%f_eval + 1
       if (fstatus > 0) goto 4020
       call eval_J(jstatus, n, m, X, w%J, params)
       info%g_eval = info%g_eval + 1
       if (jstatus > 0) goto 4010

       if (control%relative_tr_radius == 1) then 
          ! first, let's get diag(J^TJ)
          Jmax = 0.0
          do i = 1, n
             ! note:: assumes column-storage of J
             JtJdiag = norm2( w%J( (i-1)*m + 1 : i*m ) )
             if (JtJdiag > Jmax) Jmax = JtJdiag
          end do
          w%Delta = control%initial_radius_scale * (Jmax**2)
          if (control%print_level .ge. 3) write(control%out,3110) w%Delta
       else
          w%Delta = control%initial_radius
       end if
              
       if ( calculate_svd_J ) then
          call get_svd_J(n,m,w%J,s1,sn,control,svdstatus)
          if ((svdstatus .ne. 0).and.(control%print_level .ge. 3)) then 
             write( control%out, 3000 ) svdstatus
          end if
       end if

       normF = norm2(w%f)
       w%normF0 = normF

       !    g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g)
       w%g = -w%g
       normJF = norm2(w%g)
       w%normJF0 = normJF
       if (control%model == 8 .or. control%model == 9) w%normJFold = normJF

       if (control%model == 9) then
          ! make this relative....
          w%hybrid_tol = control%hybrid_tol * ( normJF/(0.5*(normF**2)) )
       end if
       
       ! save some data 
       info%obj = 0.5 * ( normF**2 )
       info%norm_g = normJF

       if (control%output_progress_vectors) then
          w%resvec(1) = info%obj
          w%gradvec(1) = info%norm_g
       end if
       
       select case (control%model)
       case (1) ! first-order
          w%hf(1:n**2) = zero
       case (2) ! second order
          if ( control%exact_second_derivatives ) then
             call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
             info%h_eval = info%h_eval + 1
             if (hfstatus > 0) goto 4030
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = zero
          end if
       case (3) ! barely second order (identity hessian)
          w%hf(1:n**2) = zero
          w%hf((/ ( (i-1)*n + i, i = 1,n ) /)) = one
       case (7) ! hybrid
          ! first call, so always first-derivatives only
          w%hf(1:n**2) = zero
       case (8) ! hybrid II
          ! always second order for first call...
                    if ( control%exact_second_derivatives ) then
             call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
             info%h_eval = info%h_eval + 1
             if (hfstatus > 0) goto 4030
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = zero
          end if
          w%use_second_derivatives = .true.
       case (9) ! hybrid (MNT)
          ! use first-order method initially
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
       case default
          goto 4040 ! unsupported model -- return to user
       end select
       
       
    end if


    w%iter = w%iter + 1
    if ( control%print_level >= 3 )  write( control%out , 3030 ) w%iter
    info%iter = w%iter
    
    rho  = -one ! intialize rho as a negative value
    success = .false.
    no_reductions = 0

    do while (.not. success) ! loop until successful
       no_reductions = no_reductions + 1
       if (no_reductions > max_tr_decrease+1) goto 4050

       !+++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the step                        !
       !    d                                      !   
       ! that the model thinks we should take next !
       !+++++++++++++++++++++++++++++++++++++++++++!
       call calculate_step(w%J,w%f,w%hf,w%g,n,m,w%Delta,w%d,control,info,& 
            w%calculate_step_ws)
       if (info%status .ne. 0) goto 4000
       
       !++++++++++++++++++!
       ! Accept the step? !
       !++++++++++++++++++!
       w%Xnew = X + w%d
       call eval_F(fstatus, n, m, w%Xnew, w%fnew, params)
       info%f_eval = info%f_eval + 1
       if (fstatus > 0) goto 4020
       normFnew = norm2(w%fnew)
       
       !++++++++++++++++++++++++++++!
       ! Get the value of the model !
       !      md :=   m_k(d)        !
       ! evaluated at the new step  !
       !++++++++++++++++++++++++++++!
       call evaluate_model(w%f,w%J,w%hf,w%d,md,m,n,control,w%evaluate_model_ws)
       
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the quantity                                   ! 
       !   rho = 0.5||f||^2 - 0.5||fnew||^2 =   actual_reduction  !
       !         --------------------------   ------------------- !
       !             m_k(0)  - m_k(d)         predicted_reduction !
       !                                                          !
       ! if model is good, rho should be close to one             !
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       call calculate_rho(normF,normFnew,md,rho)
       if (rho > control%eta_successful) success = .true.
       
       !++++++++++++++++++++++!
       ! Update the TR radius !
       !++++++++++++++++++++++!
       call update_trust_region_radius(rho,control,w%Delta)
       
       !+++++++++++++++++++++++!
       ! Add tests for model=8 !
       !+++++++++++++++++++++++!
       model8_success: if ( success .and. (control%model == 8) ) then

          if (.not. control%exact_second_derivatives) then 
             ! first, let's save some old values...
             w%g_old = w%g
             ! g_mixed = J_k^T r_{k+1}
             call mult_Jt(w%J,n,m,w%fnew,w%g_mixed)
          end if

          ! evaluate J and hf at the new point
          call eval_J(jstatus, n, m, w%Xnew, w%J, params)
          info%g_eval = info%g_eval + 1
          if (jstatus > 0) goto 4010
          if ( calculate_svd_J ) then
             call get_svd_J(n,m,w%J,s1,sn,control,svdstatus)
             if ((svdstatus .ne. 0).and.(control%print_level > 3)) then 
                write( control%out, 3000 ) svdstatus
             end if
          end if
          

          ! g = -J^Tf
          call mult_Jt(w%J,n,m,w%fnew,w%g)
          w%g = -w%g
       
          normF = normFnew
          normJF = norm2(w%g)
          
          decrease_grad: if (normJF > w%normJFold) then
             ! no reduction in residual...
             which_time_round: if (w%use_second_derivatives) then
                w%use_second_derivatives = .false.
                w%hf(1:n**2) = zero
                w%normJF_Newton = normJF
                success = .false.
                ! copy current vectors....
                ! (maybe streamline for production?!)
                w%fNewton(:) = w%fnew(:)
                w%JNewton(:) = w%J(:)
                w%XNewton(:) = w%Xnew(:)
             else  
                ! Gauss-Newton gave no benefit either....
                w%use_second_derivatives = .true.
                Newton_better: if ( w%normJF_Newton < normJF ) then
                   ! Newton values were better...replace
                   w%fnew(:) = w%fNewton(:)
                   w%J(:) = w%JNewton(:)
                   w%Xnew(:) = w%XNewton(:)
                   w%normJFold = w%normJF_Newton
                   normJF = w%normJF_Newton
                else
                   w%normJFold = normJF
                end if Newton_better
                success = .true.
             end if which_time_round
          else 
             w%normJFold = normJF
          end if decrease_grad
                 
       end if model8_success

       if (.not. success) then
          ! finally, check d makes progress
          if ( norm2(w%d) < epsmch * norm2(w%Xnew) ) goto 4060
       end if
    end do
    ! if we reach here, a successful step has been found
    
    ! update X and f
    X(:) = w%Xnew(:)
    w%f(:) = w%fnew(:)
    
    if ( control%model .ne. 8 ) then 
       
       if (.not. control%exact_second_derivatives) then 
          ! first, let's save some old values...
          ! g_old = -J_k^T r_k
          w%g_old = w%g
          ! g_mixed = -J_k^T r_{k+1}
          call mult_Jt(w%J,n,m,w%fnew,w%g_mixed)
          w%g_mixed = -w%g_mixed
       end if

       ! evaluate J and hf at the new point
       call eval_J(jstatus, n, m, X, w%J, params)
       info%g_eval = info%g_eval + 1
       if (jstatus > 0) goto 4010
       if ( calculate_svd_J ) then
          call get_svd_J(n,m,w%J,s1,sn,control,svdstatus)
          if ((svdstatus .ne. 0).and.(control%print_level > 3)) then 
             write( control%out, 3000 ) svdstatus
          end if
       end if
       
       ! g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g)
       w%g = -w%g

       if ( control%model == 9) w%normJFold = normJF

       normF = normFnew
       normJF = norm2(w%g)
       
    end if

    ! setup the vectors needed if second derivatives are not available
    if (.not. control%exact_second_derivatives) then 
       
       w%y       = w%g_old   - w%g
       w%y_sharp = w%g_mixed - w%g

    end if

    select case (control%model) ! only update hessians than change..
    case (1) ! first-order
       continue
    case (2) ! second order
       call apply_second_order_info(hfstatus,n,m, & 
            X,w%f,w%hf,eval_Hf, &
            w%d, w%y, w%y_sharp,  &
            params,control,info)
!       call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
!       info%h_eval = info%h_eval + 1
       if (hfstatus > 0) goto 4030
    case (3) ! barely second order (identity hessian)
       continue
    case (7) ! hybrid method
       if ( w%use_second_derivatives ) then 
          ! hybrid switch turned on....
          
          call apply_second_order_info(hfstatus,n,m, &
               X,w%f,w%hf,eval_Hf, &
               w%d, w%y, w%y_sharp,  &
               params,control,info)
!          call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
!          info%h_eval = info%h_eval + 1
          if (hfstatus > 0) goto 4030
       end if
    case (8)       
       call apply_second_order_info(hfstatus,n,m,&
            X,w%f,w%hf,eval_Hf, &
            w%d, w%y, w%y_sharp,  &
            params,control,info)
!       call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
!       info%h_eval = info%h_eval + 1
       if (hfstatus > 0) goto 4030
    case (9)
       ! First, check if we need to switch methods
       if (w%use_second_derivatives) then 
          if (normJf > w%normJFold) then 
             ! switch to Gauss-Newton             
             if (control%print_level .ge. 3) write(control%out,3120) 
             w%use_second_derivatives = .false.
          end if
       else
          FunctionValue = 0.5 * (normF**2)
          if ( normJf/FunctionValue < w%hybrid_tol ) then 
             w%hybrid_count = w%hybrid_count + 1
             if (w%hybrid_count == control%hybrid_switch_its) then
                ! use (Quasi-)Newton
                if (control%print_level .ge. 3) write(control%out,3130) 
                w%use_second_derivatives = .true.
                w%hybrid_count = 0
             end if
          else 
             w%hybrid_count = 0
          end if
       end if

       if (w%use_second_derivatives) then 
          call apply_second_order_info(hfstatus,n,m, &
               X,w%f,w%hf,eval_Hf, &
               w%d, w%y, w%y_sharp,  &
               params,control,info)
!          call eval_HF(hfstatus, n, m, X, w%f, w%hf, params)
 !         info%h_eval = info%h_eval + 1
          if (hfstatus > 0) goto 4030
       else 
          w%hf(1:n**2) = zero
       end if
    end select

    ! update the stats 
    info%obj = 0.5*(normF**2)
    info%norm_g = normJF
    if (control%output_progress_vectors) then
       w%resvec (w%iter + 1) = info%obj
       w%gradvec(w%iter + 1) = info%norm_g
    end if
    
    if ( (control%model == 7) .and. (normJF < control%hybrid_switch) ) then
       w%use_second_derivatives = .true.
    end if
    
    if (control%print_level >=3) write(control%out,3010) info%obj
    if (control%print_level >=3) write(control%out,3060) normJF/normF

    !++++++++++++++++++!
    ! Test convergence !
    !++++++++++++++++++!
    call test_convergence(normF,normJF,w%normF0,w%normJF0,control,info)
    if (info%convergence_normf == 1) goto 5000 ! <----converged!!
    if (info%convergence_normg == 1) goto 5010 ! <----converged!!

    if (control%print_level > 2 ) write(control%out,3100) rho

! Non-executable statements

! print level > 0

!1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

! print level > 1

! print level > 2
3000 FORMAT(/,'* Running RAL_NLLS *')
3010 FORMAT('0.5 ||f||^2 = ',ES12.4)
3030 FORMAT('== Starting iteration ',i0,' ==')
3060 FORMAT('||J''f||/||f|| = ',ES12.4)

!3090 FORMAT('Step was unsuccessful -- rho =', ES12.4)
3100 FORMAT('Step was successful -- rho =', ES12.4)
3110 FORMAT('Initial trust region radius taken as ', ES12.4)
3120 FORMAT('** Switching to Gauss-Newton **')
3130 FORMAT('** Switching to (Quasi-)Newton **')
! error returns
4000 continue
    ! generic end of algorithm
        ! all (final) exits should pass through here...
    if (control%output_progress_vectors) then
       if(.not. allocated(info%resvec)) then 
          allocate(info%resvec(w%iter + 1), stat = astatus)
          if (astatus > 0) then
             info%status = -9999
             return
          end if
          info%resvec(1:w%iter + 1) = w%resvec(1:w%iter + 1)
       end if
       if(.not. allocated(info%gradvec)) then 
          allocate(info%gradvec(w%iter + 1), stat = astatus)
          if (astatus > 0) then
             info%status = -9999
             return
          end if
          info%gradvec(1:w%iter + 1) = w%gradvec(1:w%iter + 1)
       end if
    end if

    return

4010 continue
    ! Error in eval_J
    if (control%print_level > 0) then
       write(control%error,'(a,i0)') 'Error code from eval_J, status = ', jstatus
    end if
    info%status = -2
    goto 4000

4020 continue
    ! Error in eval_J
    if (control%print_level > 0) then
       write(control%error,'(a,i0)') 'Error code from eval_F, status = ', fstatus
    end if
    info%status = -2
    goto 4000

4030 continue
    ! Error in eval_HF
    if (control%print_level > 0) then
       write(control%error,'(a,i0)') 'Error code from eval_HF, status = ', hfstatus
    end if
    info%status = -2
    goto 4000

4040 continue 
    ! unsupported choice of model
    if (control%print_level > 0) then
       write(control%error,'(a,i0,a)') 'Error: the choice of control%model = ', &
            control%model, ' is not supported'
    end if
    info%status = -3
   goto 4000

4050 continue 
    ! max tr reductions exceeded
    if (control%print_level > 0) then
       write(control%error,'(a)') 'Error: maximum tr reductions reached'
    end if
    info%status = -500
    goto 4000

4060 continue 
    ! x makes no progress
    if (control%print_level > 0) then
       write(control%error,'(a)') 'No further progress in X'
    end if
    info%status = -700
    goto 4000

4070 continue
    ! n > m on entry
    if (control%print_level > 0) then
       write(control%error,'(a)') ''
    end if
    info%status = -800
    goto 4000

! convergence 
5000 continue
    ! convegence test satisfied
    if (control%print_level > 2) then
       write(control%out,'(a,i0)') 'RAL_NLLS converged (on ||f|| test) at iteration ', &
            w%iter
    end if
    goto 4000

5010 continue
    if (control%print_level > 2) then
       write(control%out,'(a,i0)') 'RAL_NLLS converged (on gradient test) at iteration ', &
            w%iter
    end if
    goto 4000


!  End of subroutine RAL_NLLS_iterate

  end subroutine ral_nlls_iterate


  SUBROUTINE calculate_step(J,f,hf,g,n,m,Delta,d,control,info,w)

! -------------------------------------------------------
! calculate_step, find the next step in the optimization
! -------------------------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: control
     TYPE( NLLS_inform_type ), INTENT( INOUT ) :: info
     TYPE( calculate_step_work ) :: w
     select case (control%nlls_method)
        
     case (1) ! Powell's dogleg
        call dogleg(J,f,hf,g,n,m,Delta,d,control,info,w%dogleg_ws)
     case (2) ! The AINT method
        call AINT_TR(J,f,hf,n,m,Delta,d,control,info,w%AINT_tr_ws)
     case (3) ! More-Sorensen
        call more_sorensen(J,f,hf,n,m,Delta,d,control,info,w%more_sorensen_ws)
     case (4) ! Galahad
        call solve_dtrs(J,f,hf,n,m,Delta,d,w%solve_dtrs_ws,control,info)
     case default
        
        if ( control%print_level > 0 ) then
           write(control%error,'(a)') 'Error: unknown value of control%nlls_method'
           write(control%error,'(a,i0)') 'control%nlls_method = ', control%nlls_method
           info%status = -110 ! fix me
        end if

     end select

   END SUBROUTINE calculate_step


   SUBROUTINE dogleg(J,f,hf,g,n,m,Delta,d,control,info,w)
! -----------------------------------------
! dogleg, implement Powell's dogleg method
! -----------------------------------------

     REAL(wp), intent(in) :: J(:), hf(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: control
     TYPE( NLLS_inform_type ), INTENT( INOUT ) :: info
     TYPE( dogleg_work ) :: w
     
     real(wp) :: alpha, beta
     integer :: slls_status, fb_status

     !     Jg = J * g
     call mult_J(J,n,m,g,w%Jg)

     alpha = norm2(g)**2 / norm2( w%Jg )**2
       
     w%d_sd = alpha * g;

     ! Solve the linear problem...
     select case (control%model)
     case (1)
        ! linear model...
        call solve_LLS(J,f,n,m,control%lls_solver,w%d_gn,slls_status,w%solve_LLS_ws)
        if ( slls_status .ne. 0 ) goto 1000
     case default
        if (control%print_level> 0) then
           write(control%error,'(a)') 'Error: model not supported in dogleg'
        end if
        info%status = -3
        return
     end select
     
     if (norm2(w%d_gn) <= Delta) then
        d = w%d_gn
     else if (norm2( alpha * w%d_sd ) >= Delta) then
        d = (Delta / norm2(w%d_sd) ) * w%d_sd
     else
        w%ghat = w%d_gn - alpha * w%d_sd
        call findbeta(w%d_sd,w%ghat,alpha,Delta,beta,fb_status)
        if ( fb_status .ne. 0 ) goto 1010
        d = alpha * w%d_sd + beta * w%ghat
     end if
     
     return
     
1000 continue 
     ! bad error return from solve_LLS
     if ( control%print_level > 0 ) then
        write(control%out,'(a,a)') 'Unexpected error in solving a linear least squares ', &
                                   'problem in dogleg'
        write(control%out,'(a,i0)') 'dposv returned info = ', slls_status
     end if
     info%status = -700
     return

1010 continue
          if ( control%print_level > 0 ) then
        write(control%out,'(a,a)') 'Unexpected error in finding beta ', &
                                   'in dogleg'
        write(control%out,'(a,i0)') 'dogleg returned info = ', fb_status
     end if
     info%status = -701
     return

     
   END SUBROUTINE dogleg
     
   SUBROUTINE AINT_tr(J,f,hf,n,m,Delta,d,control,info,w)
     ! -----------------------------------------
     ! AINT_tr
     ! Solve the trust-region subproblem using 
     ! the method of ADACHI, IWATA, NAKATSUKASA and TAKEDA
     ! -----------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: control
     TYPE( NLLS_inform_type ), INTENT( INOUT ) :: info
     type( AINT_tr_work ) :: w
        
     integer :: solve_status, find_status
     integer :: keep_p0, i, eig_info, size_hard(2)
     real(wp) :: obj_p0, obj_p1
     REAL(wp) :: norm_p0, tau, lam, eta
     REAL(wp), allocatable :: y_hardcase(:,:)
     
     ! todo..
     ! seems wasteful to have a copy of A and B in M0 and M1
     ! use a pointer?

     keep_p0 = 0
     tau = 1e-4
     obj_p0 = HUGE(wp)

     ! The code finds 
     !  min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p||_B \leq Delta
     !
     ! set A and v for the model being considered here...

     call matmult_inner(J,n,m,w%A)
     ! add any second order information...
     w%A = w%A + reshape(hf,(/n,n/))
     call mult_Jt(J,n,m,f,w%v)
        
     ! Set B to I by hand  
     ! todo: make this an option
     w%B = 0
     do i = 1,n
        w%B(i,i) = 1.0
     end do
     
     select case (control%model)
     case (1,3)
        call solve_spd(w%A,-w%v,w%LtL,w%p0,n,solve_status)!,w%solve_spd_ws)
        if (solve_status .ne. 0) goto 1010
     case default
       call solve_general(w%A,-w%v,w%p0,n,solve_status,w%solve_general_ws)
       if (solve_status .ne. 0) goto 1020
     end select
          
     call matrix_norm(w%p0,w%B,norm_p0)
     
     if (norm_p0 < Delta) then
        keep_p0 = 1;
        ! get obj_p0 : the value of the model at p0
        call evaluate_model(f,J,hf,w%p0,obj_p0,m,n,control,w%evaluate_model_ws)
     end if

     w%M0(1:n,1:n) = -w%B
     w%M0(n+1:2*n,1:n) = w%A
     w%M0(1:n,n+1:2*n) = w%A
     call outer_product(w%v,n,w%gtg) ! gtg = Jtg * Jtg^T
     w%M0(n+1:2*n,n+1:2*n) = (-1.0 / Delta**2) * w%gtg

     w%M1 = 0.0
     w%M1(n+1:2*n,1:n) = -w%B
     w%M1(1:n,n+1:2*n) = -w%B
     
     call max_eig(w%M0,w%M1,2*n,lam, w%y, eig_info, y_hardcase,  control, w%max_eig_ws)
     if ( eig_info > 0 ) goto 1030

     if (norm2(w%y(1:n)) < tau) then
        ! Hard case
        ! overwrite H onto M0, and the outer prod onto M1...
        size_hard = shape(y_hardcase)
        call matmult_outer( matmul(w%B,y_hardcase), size_hard(2), n, w%M1(1:n,1:n))
        w%M0(1:n,1:n) = w%A(:,:) + lam*w%B(:,:) + w%M1(1:n,1:n)
        ! solve Hq + g = 0 for q
        select case (control%model) 
        case (1,3)
           call solve_spd(w%M0(1:n,1:n),-w%v,w%LtL,w%q,n,solve_status)!,w%solve_spd_ws)
        case default
          call solve_general(w%M0(1:n,1:n),-w%v,w%q,n,solve_status,w%solve_general_ws)
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! (I think..) and inside...fix

        
        ! find max eta st ||q + eta v(:,1)||_B = Delta
        call findbeta(w%q,y_hardcase(:,1),one,Delta,eta,find_status)
        if ( find_status .ne. 0 ) goto 1040

        !!!!!      ^^TODO^^    !!!!!
        ! currently assumes B = I !!
        !!!!       fixme!!      !!!!
        
        w%p1(:) = w%q(:) + eta * y_hardcase(:,1)
        
     else 
        select case (control%model)
        case (1,3)
           call solve_spd(w%A + lam*w%B,-w%v,w%LtL,w%p1,n,solve_status)!,w%solve_spd_ws)
        case default
           call solve_general(w%A + lam*w%B,-w%v,w%p1,n,solve_status,w%solve_general_ws)
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! and inside...fix
     end if
     
     ! get obj_p1 : the value of the model at p1
     call evaluate_model(f,J,hf,w%p1,obj_p1,m,n,control,w%evaluate_model_ws)

     ! what gives the smallest objective: p0 or p1?
     if (obj_p0 < obj_p1) then
        d = w%p0
     else 
        d = w%p1
     end if

     return

1010 continue 
     ! bad error return from solve_spd
     if ( control%print_level >= 0 ) then 
        write(control%error,'(a)') 'Error in solving a linear system in AINT_TR'
        write(control%error,'(a,i0)') 'dposv returned info = ', solve_status
     end if
     info%status = -4
     return
     
1020 continue
     ! bad error return from solve_general
     if ( control%print_level >= 0 ) then 
        write(control%error,'(a)') 'Error in solving a linear system in AINT_TR'
        write(control%error,'(a,i0)') 'dgexv returned info = ', solve_status
     end if
     info%status = -4
     return
     
1030 continue
     ! bad error return from max_eig
     if ( control%print_level >= 0 ) then 
        write(control%error,'(a)') 'Error in the eigenvalue computation of AINT_TR'
        write(control%error,'(a,i0)') 'dggev returned info = ', eig_info
     end if
     info%status = -4
     return

1040 continue
     ! no valid beta found
     if ( control%print_level >= 0 ) then 
        write(control%error,'(a)') 'No valid beta found'
     end if
     info%status = -4
     return

   END SUBROUTINE AINT_tr

   subroutine more_sorensen(J,f,hf,n,m,Delta,d,control,info,w)
     ! -----------------------------------------
     ! more_sorensen
     ! Solve the trust-region subproblem using 
     ! the method of More and Sorensen
     !
     ! Using the implementation as in Algorithm 7.3.6
     ! of Trust Region Methods
     ! 
     ! main output :: d, the soln to the TR subproblem
     ! -----------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     TYPE( NLLS_control_type ), INTENT( IN ) :: control
     TYPE( NLLS_inform_type ), INTENT( INOUT ) :: info
     type( more_sorensen_work ) :: w

     ! parameters...make these control?
     real(wp) :: nd, nq

     real(wp) :: sigma, alpha
     integer :: fb_status, mineig_status
     integer :: test_pd, i, no_shifts
     
     ! The code finds 
     !  d = arg min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p|| \leq Delta
     !
     ! set A and v for the model being considered here...
     
     ! Set A = J^T J
     call matmult_inner(J,n,m,w%A)
     ! add any second order information...
     ! so A = J^T J + HF
     w%A = w%A + reshape(hf,(/n,n/))
     ! now form v = J^T f 
     call mult_Jt(J,n,m,f,w%v)
          
     ! d = -A\v
     call solve_spd(w%A,-w%v,w%LtL,d,n,test_pd)!,w%solve_spd_ws)
     if (test_pd .eq. 0) then
        ! A is symmetric positive definite....
        sigma = zero
     else
        call min_eig_symm(w%A,n,sigma,w%y1,control,mineig_status,w%min_eig_symm_ws) 
        if (mineig_status .ne. 0) goto 1060 
        sigma = -(sigma - control%more_sorensen_shift)
        no_shifts = 1
100     call shift_matrix(w%A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-w%v,w%LtL,d,n,test_pd)
        if ( test_pd .ne. 0 ) then
           no_shifts = no_shifts + 1
           if ( no_shifts == 10 ) goto 3000
           sigma =  sigma + (10**no_shifts) * control%more_sorensen_shift
           if (control%print_level >=3) write(control%out,2000) sigma
           goto 100 
        end if
     end if
     
     nd = norm2(d)
     if (nd .le. Delta) then
        ! we're within the tr radius from the start!
        if ( abs(sigma) < control%more_sorensen_tiny ) then
           ! we're good....exit
           goto 1050
        else if ( abs( nd - Delta ) < control%more_sorensen_tiny ) then
           ! also good...exit
           goto 1050              
        end if
        call findbeta(d,w%y1,one,Delta,alpha,fb_status)
        if (fb_status .ne. 0 ) goto 1070  !! todo! change this error code....
        d = d + alpha * w%y1
        ! also good....exit
        goto 1050
     end if

     ! now, we're not in the trust region initally, so iterate....
     do i = 1, control%more_sorensen_maxits
        if ( abs(nd - Delta) <= control%more_sorensen_tol * Delta) then
           goto 1020 ! converged!
        end if
        
        w%q = d ! w%q = R'\d
        CALL DTRSM( 'Left', 'Lower', 'No Transpose', 'Non-unit', n, & 
             1, one, w%LtL, n, w%q, n )
        
        nq = norm2(w%q)
        
        sigma = sigma + ( (nd/nq)**2 )* ( (nd - Delta) / Delta )
        
        call shift_matrix(w%A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-w%v,w%LtL,d,n,test_pd)
        if (test_pd .ne. 0)  goto 2010 ! shouldn't happen...
        
        nd = norm2(d)

     end do
     
     goto 1040
     
1020 continue
     ! Converged!
     if ( control%print_level >= 3 ) then
        write(control%error,'(a,i0)') 'More-Sorensen converged at iteration ', i
     end if
     return
     
1040 continue
     ! maxits reached, not converged
     if ( control%print_level > 0 ) then
        write(control%error,'(a)') 'Maximum iterations reached in More-Sorensen'
        write(control%error,'(a)') 'without convergence'
     end if
     info%status = -100 ! fix me
     return

1050 continue
     if ( control%print_level >= 3 ) then
        write(control%error,'(a)') 'More-Sorensen: first point within trust region'
     end if
     return

1060 continue
     if ( control%print_level > 0 ) then
        write(control%error,'(a)') 'More-Sorensen: error from lapack routine dsyev(x)'
        write(control%error,'(a,i0)') 'info = ', mineig_status
     end if
     info%status = -333

     return

1070 continue
     if ( control%print_level >= 3 ) then
        write(control%error,'(a)') 'M-S: Unable to find alpha s.t. ||s + alpha v|| = Delta'
     end if
     info%status = -200

     return
     
2000 format('Non-spd system in more_sorensen. Increasing sigma to ',es12.4)

3000 continue 
     ! bad error return from solve_spd
     if ( control%print_level > 0 ) then
        write(control%out,'(a)') 'Unexpected error in solving a linear system in More_sorensen'
        write(control%out,'(a,i0)') 'dposv returned info = ', test_pd
     end if
     info%status = -500
     return
     
2010 continue 
     ! bad error return from solve_spd
     if ( control%print_level > 0 ) then
        write(control%out,'(a,a)') 'Unexpected error in solving a linear system ', &
                                   'in More_sorensen loop'
        write(control%out,'(a,i0)') 'dposv returned info = ', test_pd
     end if
     info%status = -600
     return
     
     
   end subroutine more_sorensen
   
   subroutine solve_dtrs(J,f,hf,n,m,Delta,d,w,control,info)

     !---------------------------------------------
     ! solve_dtrs
     ! Solve the trust-region subproblem using
     ! the DTRS method from Galahad
     ! 
     ! This method needs H to be diagonal, so we need to 
     ! pre-process
     !
     ! main output :: d, the soln to the TR subproblem
     !--------------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     type( solve_dtrs_work ) :: w
     TYPE( NLLS_control_type ), INTENT( IN ) :: control
     TYPE( NLLS_inform_type ), INTENT( INOUT ) :: info

     integer :: eig_status
     TYPE ( DTRS_CONTROL_TYPE ) :: dtrs_control
     TYPE ( DTRS_inform_type )  :: dtrs_inform
     
     ! The code finds 
     !  d = arg min_p   w^T p + 0.5 * p^T D p
     !       s.t. ||p|| \leq Delta
     !
     ! where D is diagonal
     !
     ! our probem in naturally in the form
     ! 
     ! d = arg min_p   v^T p + 0.5 * p^T H p
     !       s.t. ||p|| \leq Delta
     !
     ! first, find the matrix H and vector v
     ! Set A = J^T J
     call matmult_inner(J,n,m,w%A)
     ! add any second order information...
     ! so A = J^T J + HF
     w%A = w%A + reshape(hf,(/n,n/))

     ! now form v = J^T f 
     call mult_Jt(J,n,m,f,w%v)

     ! Now that we have the unprocessed matrices, we need to get an 
     ! eigendecomposition to make A diagonal
     !
     call all_eig_symm(w%A,n,w%ew,w%ev,w%all_eig_symm_ws,eig_status)
     if (eig_status .ne. 0) goto 1000

     ! We can now change variables, setting y = Vp, getting
     ! Vd = arg min_(Vx) v^T p + 0.5 * (Vp)^T D (Vp)
     !       s.t.  ||x|| \leq Delta
     ! <=>
     ! Vd = arg min_(Vx) V^Tv^T (Vp) + 0.5 * (Vp)^T D (Vp)
     !       s.t.  ||x|| \leq Delta
     ! <=>

     ! we need to get the transformed vector v
     call mult_Jt(w%ev,n,n,w%v,w%v_trans)

     ! we've now got the vectors we need, pass to dtrs_solve
     call dtrs_initialize( dtrs_control, dtrs_inform ) 

     call dtrs_solve(n, Delta, zero, w%v_trans, w%ew, w%d_trans, dtrs_control, dtrs_inform )
     if ( dtrs_inform%status .ne. 0) goto 1010

     ! and return the un-transformed vector
     call mult_J(w%ev,n,n,w%d_trans,d)

     return
     
1000 continue
     if ( control%print_level > 0 ) then
        write(control%error,'(a)') 'solve_dtrs: error from lapack routine dsyev(x)'
        write(control%error,'(a,i0)') 'info = ', eig_status
     end if
     info%status = -333
     return

1010 continue
     if ( control%print_level > 0 ) then
        write(control%error,'(a)') 'solve_dtrs: error from GALAHED routine DTRS'
        write(control%error,'(a,i0)') 'info = ', dtrs_inform%status
     end if
     info%status = -777
     return

   end subroutine solve_dtrs


   SUBROUTINE solve_LLS(J,f,n,m,method,d_gn,status,w)
       
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
       type( solve_LLS_work ) :: w
       
       lda = m
       ldb = max(m,n)
       w%temp(1:m) = f(1:m)
       lwork = size(w%work)
       
       w%Jlls(:) = J(:)

       call dgels(trans, m, n, nrhs, w%Jlls, lda, w%temp, ldb, w%work, lwork, status)
       
       d_gn = -w%temp(1:n)
              
     END SUBROUTINE solve_LLS
     
     SUBROUTINE findbeta(d_sd,ghat,alpha,Delta,beta,status)

!  -----------------------------------------------------------------
!  findbeta, a subroutine to find the optimal beta such that 
!   || d || = Delta, where d = alpha * d_sd + beta * ghat
!  -----------------------------------------------------------------

     real(wp), dimension(:), intent(in) :: d_sd, ghat
     real(wp), intent(in) :: alpha, Delta
     real(wp), intent(out) :: beta
     integer, intent(out) :: status
     
     real(wp) :: a, b, c, discriminant
     
     status = 0

     a = norm2(ghat)**2
     b = 2.0 * alpha * dot_product( ghat, d_sd)
     c = ( alpha * norm2( d_sd ) )**2 - Delta**2
     
     discriminant = b**2 - 4 * a * c
     if ( discriminant < 0) then
        status = 1
        return
     else
        beta = (-b + sqrt(discriminant)) / (2.0 * a)
     end if

     END SUBROUTINE findbeta

     
     subroutine evaluate_model(f,J,hf,d,md,m,n,control,w)
! --------------------------------------------------
! Input:
! f = f(x_k), J = J(x_k), 
! hf = \sum_{i=1}^m f_i(x_k) \nabla^2 f_i(x_k) (or an approx)
!
! We have a model 
!      m_k(d) = 0.5 f^T f  + d^T J f + 0.5 d^T (J^T J + HF) d
!
! This subroutine evaluates the model at the point d 
! This value is returned as the scalar
!       md :=m_k(d)
! --------------------------------------------------       

       real(wp), intent(in) :: f(:) ! f(x_k)
       real(wp), intent(in) :: d(:) ! direction in which we move
       real(wp), intent(in) :: J(:) ! J(x_k) (by columns)
       real(wp), intent(in) :: hf(:)! (approx to) \sum_{i=1}^m f_i(x_k) \nabla^2 f_i(x_k)
       integer, intent(in) :: m,n
       real(wp), intent(out) :: md  ! m_k(d)
       TYPE( NLLS_control_type ), INTENT( IN ) :: control
       type( evaluate_model_work ) :: w
       
       !Jd = J*d
       call mult_J(J,n,m,d,w%Jd)
       
       ! First, get the base 
       ! 0.5 (f^T f + f^T J d + d^T' J ^T J d )
       md = 0.5 * norm2(f + w%Jd)**2
       select case (control%model)
       case (1) ! first-order (no Hessian)
          ! nothing to do here...
          continue
       case (3) ! barely second-order (identity Hessian)
          ! H = J^T J + I
          md = md + 0.5 * dot_product(d,d)
       case default
          ! these have a dynamic H -- recalculate
          ! H = J^T J + HF, HF is (an approx?) to the Hessian
          call mult_J(hf,n,n,d,w%Hd)
          md = md + 0.5 * dot_product(d,w%Hd)
       end select

     end subroutine evaluate_model
     
     subroutine calculate_rho(normf,normfnew,md,rho)
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the quantity                                   ! 
       !   rho = 0.5||f||^2 - 0.5||fnew||^2 =   actual_reduction  !
       !         --------------------------   ------------------- !
       !             m_k(0)  - m_k(d)         predicted_reduction !
       !                                                          !
       ! if model is good, rho should be close to one             !
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       real(wp), intent(in)  :: normf    ! ||f(x_k)|| at 
       real(wp), intent(in)  :: normfnew ! ||f(x_k + d)||
       real(wp), intent(in)  :: md       !    m_k(d)
       real(wp), intent(out) :: rho      ! act_red / pred_red (close to 1 == good)

       real(wp) :: actual_reduction, predicted_reduction
       
       actual_reduction = ( 0.5 * normf**2 ) - ( 0.5 * normfnew**2 )
       predicted_reduction = ( ( 0.5 * normf**2 ) - md )
       
       if ( abs(actual_reduction) < 10*epsmch ) then 
          rho = one
       else if (abs( predicted_reduction ) < 10 * epsmch ) then 
          rho = one
       else
          rho = actual_reduction / predicted_reduction
       end if


     end subroutine calculate_rho

     subroutine apply_second_order_info(status,n,m,&
          X,f,hf,eval_Hf,&
          d, y, y_sharp, & 
          params,control,info)
       integer, intent(out) :: status
       integer, intent(in)  :: n, m 
       real(wp), intent(in) :: X(:), f(:)
       real(wp), intent(inout) :: hf(:)
       real(wp), intent(in) :: d(:), y(:), y_sharp(:)
       procedure( eval_hf_type ) :: eval_Hf
       class( params_base_type ) :: params
       type( NLLS_control_type ), intent(in) :: control
       type( NLLS_inform_type ), intent(inout) :: info
       
       real(wp), allocatable :: Sks(:), ysharpSks(:)

       real(wp) :: yts, alpha
       integer :: i,j

       if (control%exact_second_derivatives) then
          call eval_HF(status, n, m, X, f, hf, params)
          info%h_eval = info%h_eval + 1
       else

          yts = dot_product(d,y)
          allocate(Sks(n))
          call mult_J(hf,n,n,d,Sks) ! hfs = S_k * d

          allocate(ysharpSks(n))
          ysharpSks = y_sharp - Sks
          
          ! now, let's scale hd (Nocedal and Wright, Section 10.2)
          alpha = abs(dot_product(d,y_sharp))/abs(dot_product(d,Sks))
          alpha = min(one,alpha)
          hf(:)  = alpha * hf(:)

          ! update S_k (again, as in N&W, Section 10.2)

          ! hf = hf + (1/yts) (y# - Sk d)^T y:
          alpha = 1/yts
          call dGER(n,n,alpha,ysharpSks,1,y,1,hf,n)
          ! hf = hf + (1/yts) y^T (y# - Sk d):
          call dGER(n,n,alpha,y,1,ysharpSks,1,hf,n)
          ! hf = hf - ((y# - Sk d)^T d)/((yts)**2)) * y y^T
          alpha = -dot_product(ysharpSks,d)/(yts**2)
          call dGER(n,n,alpha,y,1,y,1,hf,n)
                    
       end if

     end subroutine apply_second_order_info

     subroutine update_trust_region_radius(rho,control,Delta)
       
       real(wp), intent(inout) :: rho ! ratio of actual to predicted reduction
       type( NLLS_control_type ), intent(in) :: control
       real(wp), intent(inout) :: Delta ! trust region size
       
       if (rho < control%eta_successful) then
          ! unsuccessful....reduce Delta
          Delta = max( control%radius_reduce, control%radius_reduce_max) * Delta
          if (control%print_level > 2) write(control%out,3010) Delta     
       else if (rho < control%eta_very_successful) then 
          ! doing ok...retain status quo
          if (control%print_level > 2) write(control%out,3020) Delta 
       else if (rho < control%eta_too_successful ) then
          ! more than very successful -- increase delta
          Delta = min(control%maximum_radius, control%radius_increase * Delta )
          if (control%print_level > 2) write(control%out,3030) Delta
       else if (rho >= control%eta_too_successful) then
          ! too successful....accept step, but don't change Delta
          if (control%print_level > 2) write(control%out,3040) Delta 
       else
          ! just incase (NaNs and the like...)
          if (control%print_level > 2) write(control%out,3010) Delta 
          Delta = max( control%radius_reduce, control%radius_reduce_max) * Delta
          rho = -one ! set to be negative, so that the logic works....
       end if 

       return
       
! print statements

3010   FORMAT('Unsuccessful step -- decreasing Delta to', ES12.4)      
3020   FORMAT('Successful step -- Delta staying at', ES12.4)     
3030   FORMAT('Very successful step -- increasing Delta to', ES12.4)
3040   FORMAT('Step too successful -- Delta staying at', ES12.4)   
       

     end subroutine update_trust_region_radius
     
     subroutine test_convergence(normF,normJF,normF0,normJF0,control,info)
       
       real(wp), intent(in) :: normF, normJf, normF0, normJF0
       type( NLLS_control_type ), intent(in) :: control
       type( NLLS_inform_type ), intent(inout) :: info

       if ( normF <= max(control%stop_g_absolute, &
            control%stop_g_relative * normF0) ) then
          info%convergence_normf = 1
          return
       end if
       
       if ( (normJF/normF) <= max(control%stop_g_absolute, &
            control%stop_g_relative * (normJF0/normF0)) ) then
          info%convergence_normg = 1
       end if
          
       return
       
     end subroutine test_convergence
     
     subroutine mult_J(J,n,m,x,Jx)
       real(wp), intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       real(wp), intent(out) :: Jx(*)
       
       real(wp) :: alpha, beta

       Jx(1:m) = 1.0
       alpha = 1.0
       beta  = 0.0

       call dgemv('N',m,n,alpha,J,m,x,1,beta,Jx,1)
       
     end subroutine mult_J

     subroutine mult_Jt(J,n,m,x,Jtx)
       double precision, intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       double precision, intent(out) :: Jtx(*)
       
       double precision :: alpha, beta

       Jtx(1:n) = one
       alpha = one
       beta  = zero

       call dgemv('T',m,n,alpha,J,m,x,1,beta,Jtx,1)

      end subroutine mult_Jt

      subroutine solve_spd(A,b,LtL,x,n,status)
        REAL(wp), intent(in) :: A(:,:)
        REAL(wp), intent(in) :: b(:)
        REAL(wp), intent(out) :: LtL(:,:)
        REAL(wp), intent(out) :: x(:)
        integer, intent(in) :: n
        integer, intent(out) :: status

        ! A wrapper for the lapack subroutine dposv.f
        ! get workspace for the factors....
        status = 0
        LtL(1:n,1:n) = A(1:n,1:n)
        x(1:n) = b(1:n)
        call dposv('L', n, 1, LtL, n, x, n, status)
           
      end subroutine solve_spd

      subroutine solve_general(A,b,x,n,status,w)
        REAL(wp), intent(in) :: A(:,:)
        REAL(wp), intent(in) :: b(:)
        REAL(wp), intent(out) :: x(:)
        integer, intent(in) :: n
        integer, intent(out) :: status
        type( solve_general_work ) :: w
        
        ! A wrapper for the lapack subroutine dposv.f
        ! NOTE: A would be destroyed
        w%A(1:n,1:n) = A(1:n,1:n)
        x(1:n) = b(1:n)
        call dgesv( n, 1, w%A, n, w%ipiv, x, n, status)
        
      end subroutine solve_general

      subroutine matrix_norm(x,A,norm_A_x)
        REAL(wp), intent(in) :: A(:,:), x(:)
        REAL(wp), intent(out) :: norm_A_x

        ! Calculates norm_A_x = ||x||_A = sqrt(x'*A*x)

        norm_A_x = sqrt(dot_product(x,matmul(A,x)))

      end subroutine matrix_norm

      subroutine matmult_inner(J,n,m,A)
        
        integer, intent(in) :: n,m 
        real(wp), intent(in) :: J(*)
        real(wp), intent(out) :: A(n,n)
        integer :: lengthJ
        
        ! Takes an m x n matrix J and forms the 
        ! n x n matrix A given by
        ! A = J' * J
        
        lengthJ = n*m
        
        call dgemm('T','N',n, n, m, one,&
                   J, m, J, m, & 
                   zero, A, n)
        
        
      end subroutine matmult_inner

       subroutine matmult_outer(J,n,m,A)
        
        integer, intent(in) :: n,m 
        real(wp), intent(in) :: J(*)
        real(wp), intent(out) :: A(m,m)
        integer :: lengthJ

        ! Takes an m x n matrix J and forms the 
        ! m x m matrix A given by
        ! A = J * J'
        
        lengthJ = n*m
        
        call dgemm('N','T',m, m, n, one,&
                   J, m, J, m, & 
                   zero, A, m)
        
        
      end subroutine matmult_outer
      
      subroutine outer_product(x,n,xtx)
        
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: n
        real(wp), intent(out) :: xtx(:,:)

        ! Takes an n vector x and forms the 
        ! n x n matrix xtx given by
        ! xtx = x * x'

        xtx(1:n,1:n) = zero
        call dger(n, n, one, x, 1, x, 1, xtx, n)
        
      end subroutine outer_product

      subroutine all_eig_symm(A,n,ew,ev,w,status)
        ! calculate all the eigenvalues of A (symmetric)

        real(wp), intent(in) :: A(:,:)
        integer, intent(in) :: n
        real(wp), intent(out) :: ew(:), ev(:,:)
        type( all_eig_symm_work ) :: w
        integer, intent(out) :: status

        real(wp), allocatable :: work
        real(wp) :: tol
        integer :: lwork
        
        status = 0 

        ! copy the matrix A into the eigenvector array
        ev(1:n,1:n) = A(1:n,1:n)
        
        lwork = size(w%work)
        ! call dsyev --> all eigs of a symmetric matrix
        
        call dsyev('V', & ! both ew's and ev's 
             'U', & ! upper triangle of A
             n, ev, n, & ! data about A
             ew, w%work, lwork, & 
             status)
        
      end subroutine all_eig_symm

      subroutine min_eig_symm(A,n,ew,ev,control,status,w)
        ! calculate the leftmost eigenvalue of A
        
        real(wp), intent(in) :: A(:,:)
        integer, intent(in) :: n
        real(wp), intent(out) :: ew, ev(:)
        integer, intent(out) :: status
        type( NLLS_control_type ), INTENT( IN ) :: control
        type( min_eig_symm_work ) :: w

        real(wp) :: tol, dlamch
        integer :: lwork, eigsout, minindex(1)

        tol = 2*dlamch('S')!1e-15
        
        status = 0
        w%A(1:n,1:n) = A(1:n,1:n) ! copy A, as workspace for dsyev(x)
        ! note that dsyevx (but not dsyev) only destroys the lower (or upper) part of A
        ! so we could possibly reduce memory use here...leaving for 
        ! ease of understanding for now.

        lwork = size(w%work)
        if ( control%subproblem_eig_fact ) then
           ! call dsyev --> all eigs of a symmetric matrix
           call dsyev('V', & ! both ew's and ev's 
                'U', & ! upper triangle of A
                n, w%A, n, & ! data about A
                w%ew, w%work, lwork, & 
                status)
           
           minindex = minloc(w%ew)
           ew = w%ew(minindex(1))
           ev = w%A(1:n,minindex(1))
           
        else
           ! call dsyevx --> selected eigs of a symmetric matrix
           call dsyevx( 'V',& ! get both ew's and ev's
                'I',& ! just the numbered eigenvalues
                'U',& ! upper triangle of A
                n, w%A, n, & 
                1.0, 1.0, & ! not used for RANGE = 'I'
                1, 1, & ! only find the first eigenpair
                tol, & ! abstol for the eigensolver
                eigsout, & ! total number of eigs found
                ew, ev, & ! the eigenvalue and eigenvector
                n, & ! ldz (the eigenvector array)
                w%work, lwork, w%iwork, &  ! workspace
                w%ifail, & ! array containing indicies of non-converging ews
                status)

        end if
           
        ! let the calling subroutine handle the errors
        
        return
                      
      end subroutine min_eig_symm

      subroutine max_eig(A,B,n,ew,ev,status,nullevs,control,w)
        
        real(wp), intent(inout) :: A(:,:), B(:,:)
        integer, intent(in) :: n 
        real(wp), intent(out) :: ew, ev(:)
        integer, intent(out) :: status
        real(wp), intent(out), allocatable :: nullevs(:,:)
        type( NLLS_control_type ), intent(in) :: control
        type( max_eig_work ) :: w
        
        integer :: lwork, maxindex(1), no_null, halfn
        real(wp):: tau
        integer :: i 

        ! Find the max eigenvalue/vector of the generalized eigenproblem
        !     A * y = lam * B * y
        ! further, if ||y(1:n/2)|| \approx 0, find and return the 
        ! eigenvectors y(n/2+1:n) associated with this

        status = 0
        ! check that n is even (important for hard case -- see below)
        if (modulo(n,2).ne.0) goto 1010
        
        halfn = n/2
        lwork = size(w%work)
        call dggev('N', & ! No left eigenvectors
                   'V', &! Yes right eigenvectors
                   n, A, n, B, n, &
                   w%alphaR, w%alphaI, w%beta, & ! eigenvalue data
                   w%vr, n, & ! not referenced
                   w%vr, n, & ! right eigenvectors
                   w%work, lwork, status)

        ! now find the rightmost real eigenvalue
        w%vecisreal = .true.
        where ( abs(w%alphaI) > 1e-8 ) w%vecisreal = .false.
        
        w%ew_array(:) = w%alphaR(:)/w%beta(:)
        maxindex = maxloc(w%ew_array,w%vecisreal)
        if (maxindex(1) == 0) goto 1000
        
        tau = 1e-4 ! todo -- pass this through from above...
        ! note n/2 always even -- validated by test on entry
        if (norm2( w%vr(1:halfn,maxindex(1)) ) < tau) then 
           ! hard case
           ! let's find which ev's are null...
           w%nullindex = 0
           no_null = 0
           do i = 1,n
              if (norm2( w%vr(1:halfn,i)) < 1e-4 ) then
                 no_null = no_null + 1 
                 w%nullindex(no_null) = i
              end if
           end do
           allocate(nullevs(halfn,no_null))
           nullevs(:,:) = w%vr(halfn+1 : n,w%nullindex(1:no_null))
        end if
        
        ew = w%alphaR(maxindex(1))/w%beta(maxindex(1))
        ev(:) = w%vr(:,maxindex(1))

        return 

1000    continue 
        if ( control%print_level >=0 ) then
           write(control%error,'(a)') 'Error, all eigs are imaginary'
        end if
        status = 1 ! Eigs imaginary error
        
        return

1010    continue
        if (control%print_level >= 0 ) then 
           write(control%error,'(a)') 'error : non-even sized matrix sent to max eig'
        end if
        status = 2

        return

                
      end subroutine max_eig

      subroutine shift_matrix(A,sigma,AplusSigma,n)
        
        real(wp), intent(in)  :: A(:,:), sigma
        real(wp), intent(out) :: AplusSigma(:,:)
        integer, intent(in) :: n 

        integer :: i 
        ! calculate AplusSigma = A + sigma * I

        AplusSigma(:,:) = A(:,:)
        do i = 1,n
           AplusSigma(i,i) = AplusSigma(i,i) + sigma
        end do
                
      end subroutine shift_matrix

      subroutine get_svd_J(n,m,J,s1,sn,control,status)
        integer, intent(in) :: n,m 
        real(wp), intent(in) :: J(:)
        real(wp), intent(out) :: s1, sn
        type( NLLS_control_type ) :: control
        integer, intent(out) :: status

        character :: jobu(1), jobvt(1)
        real(wp), allocatable :: Jcopy(:)
        real(wp), allocatable :: S(:)
        real(wp), allocatable :: work(:)
        integer :: lwork
        
        allocate(Jcopy(n*m))
        allocate(S(n))
        Jcopy(:) = J(:)

        jobu  = 'N' ! calculate no left singular vectors
        jobvt = 'N' ! calculate no right singular vectors
        
        allocate(work(1))
        ! make a workspace query....
        call dgesvd( jobu, jobvt, n, m, Jcopy, n, S, S, 1, S, 1, & 
             work, -1, status )
        if (status .ne. 0 ) return

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))     
        
        ! call the real thing....
        call dgesvd( JOBU, JOBVT, n, m, Jcopy, n, S, S, 1, S, 1, & 
             work, lwork, status )
        if ( (status .ne. 0) .and. (control%print_level > 3) ) then 
           write(control%out,'(a,i0)') 'Error when calculating svd, dgesvd returned', &
                                        status
           s1 = -1.0
           sn = -1.0
           ! allow to continue, but warn user and return zero singular values
        else
           s1 = S(1)
           sn = S(n)
           write(control%out,'(a,es12.4,a,es12.4)') 's1 = ', s1, '    sn = ', sn
           write(control%out,'(a,es12.4)') 'k(J) = ', s1/sn
        end if

      end subroutine get_svd_J


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                       !!
!! W O R K S P A C E   S E T U P   S U B R O U T I N E S !!
!!                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_workspaces(workspace,n,m,control,status)
        
        type( NLLS_workspace ), intent(out) :: workspace
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(in) :: n,m
        integer, intent(out) :: status

        status = 0      
        
        workspace%first_call = 0

        if (.not. control%exact_second_derivatives) then
           if (.not. allocated(workspace%y)) then
              allocate(workspace%y(n), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%y_sharp)) then
              allocate(workspace%y_sharp(n), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%g_old)) then
              allocate(workspace%g_old(n), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%g_mixed)) then
              allocate(workspace%g_mixed(n), stat = status)
              if (status > 0) goto 9000
           end if

        end if

        if( control%output_progress_vectors ) then 
           if (.not. allocated(workspace%resvec)) then
              allocate(workspace%resvec(control%maxit+1), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%gradvec)) then
              allocate(workspace%gradvec(control%maxit+1), stat = status)
              if (status > 0) goto 9000
           end if
        end if

        if( control%model == 8 ) then 
           if (.not. allocated(workspace%fNewton)) then
              allocate(workspace%fNewton(m), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%JNewton)) then
              allocate(workspace%JNewton(n*m), stat = status)
              if (status > 0) goto 9000
           end if
           if (.not. allocated(workspace%XNewton)) then
              allocate(workspace%XNewton(n), stat = status)
              if (status > 0) goto 9000
           end if
        end if
                
        if( .not. allocated(workspace%J)) then
           allocate(workspace%J(n*m), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%f)) then
           allocate(workspace%f(m), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%fnew)) then 
           allocate(workspace%fnew(m), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%hf)) then
           allocate(workspace%hf(n*n), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%d)) then
           allocate(workspace%d(n), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%g)) then
           allocate(workspace%g(n), stat = status)
           if (status > 0) goto 9000
        end if
        if( .not. allocated(workspace%Xnew)) then
           allocate(workspace%Xnew(n), stat = status)
           if (status > 0) goto 9000
        end if


        select case (control%nlls_method)
        
        case (1) ! use the dogleg method
           call setup_workspace_dogleg(n,m,workspace%calculate_step_ws%dogleg_ws, & 
                control, status)
           if (status > 0) goto 9000

        case(2) ! use the AINT method
           call setup_workspace_AINT_tr(n,m,workspace%calculate_step_ws%AINT_tr_ws, & 
                control, status)
           if (status > 0) goto 9010
           
        case(3) ! More-Sorensen 
           call setup_workspace_more_sorensen(n,m,&
                workspace%calculate_step_ws%more_sorensen_ws,control,status)
           if (status > 0) goto 9000

        case (4) ! dtrs (Galahad)
           call setup_workspace_solve_dtrs(n,m, & 
                workspace%calculate_step_ws%solve_dtrs_ws, control, status)

        end select

! evaluate model in the main routine...       
        call setup_workspace_evaluate_model(n,m,workspace%evaluate_model_ws,control,status)
        if (status > 0) goto 9010

        return

! Error statements
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating local array: ',&
                'not enough memory.' 
           write(control%error,'(a,i0)') 'status = ', status
 
        end if
        return
       
9010    continue 
        
        return

      end subroutine setup_workspaces

      subroutine setup_workspace_dogleg(n,m,w,control,status)
        integer, intent(in) :: n, m 
        type( dogleg_work ), intent(out) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(inout) :: status

        allocate(w%d_sd(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%d_gn(n),stat = status)
        if (status > 0) goto 9000           
        allocate(w%ghat(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%Jg(m),stat = status)
        if (status > 0) goto 9000
        ! setup space for 
        !   solve_LLS
        call setup_workspace_solve_LLS(n,m,w%solve_LLS_ws,control,status)
        if (status > 0 ) goto 9010
        ! setup space for 
        !   evaluate_model
        call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,control,status)
        if (status > 0 ) goto 9010

        return

        ! Error statements
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''dogleg'': ',&
                'not enough memory.' 
           write(control%error,'(a,i0)') 'status = ', status
 
        end if
        
        return

9010    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a)') &
                'Called from subroutine ''dogleg'': '
        end if

        return
        

      end subroutine setup_workspace_dogleg

      subroutine setup_workspace_solve_LLS(n,m,w,control,status)
        integer, intent(in) :: n, m 
        type( solve_LLS_work ) :: w 
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(inout) :: status
        integer :: lwork
        
        allocate( w%temp(max(m,n)), stat = status)
        if (status > 0) goto 9000
        lwork = max(1, min(m,n) + max(min(m,n), 1)*4) 
        allocate( w%work(lwork), stat = status)
        if (status > 0) goto 9000
        allocate( w%Jlls(n*m), stat = status)
        if (status > 0) goto 9000
        
        return
        
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''solve_LLS'': ',&
                'not enough memory.' 
        end if
        
        return

      end subroutine setup_workspace_solve_LLS
      
      subroutine setup_workspace_evaluate_model(n,m,w,control,status)
        integer, intent(in) :: n, m        
        type( evaluate_model_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        
        allocate( w%Jd(m), stat = status )
        if (status > 0) goto 9000
        allocate( w%Hd(n), stat = status)
        if (status > 0) goto 9000

        return

9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''evaluate_model'': ',&
                'not enough memory.' 
        end if
        
        return
      end subroutine setup_workspace_evaluate_model

      subroutine setup_workspace_AINT_tr(n,m,w,control,status)
        integer, intent(in) :: n, m 
        type( AINT_tr_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        
        allocate(w%A(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%v(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%B(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%p0(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%p1(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%M0(2*n,2*n),stat = status)
        if (status > 0) goto 9000
        allocate(w%M1(2*n,2*n),stat = status)
        if (status > 0) goto 9000
        allocate(w%y(2*n),stat = status)
        if (status > 0) goto 9000
        allocate(w%gtg(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%q(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%LtL(n,n),stat = status)
        if (status > 0) goto 9000
        ! setup space for max_eig
        call setup_workspace_max_eig(n,m,w%max_eig_ws,control,status)
        if (status > 0) goto 9010
        call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,control,status)
        if (status > 0) goto 9010
        ! setup space for the solve routine
        if ((control%model .ne. 1).and.(control%model .ne. 3)) then
           call setup_workspace_solve_general(n,m,w%solve_general_ws,control,status)
           if (status > 0 ) goto 9010
        end if

        return
        
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''AINT_tr'': ',&
                'not enough memory.' 
        end if
        
        return

9010    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a)') &
                'Called from subroutine ''solve_LLS'': '
        end if
        return
        
      end subroutine setup_workspace_AINT_tr


      subroutine setup_workspace_min_eig_symm(n,m,w,control,status)
        integer, intent(in) :: n, m 
        type( min_eig_symm_work) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        
        real(wp), allocatable :: workquery(:)
        integer :: lapack_status, lwork, eigsout
        
        allocate(w%A(n,n),stat = status)
        if (status > 0) goto 9000
        
        allocate(workquery(1),stat = status)
        if (status > 0) goto 9000
        lapack_status = 0
        
        if (control%subproblem_eig_fact) then 
           allocate(w%ew(n), stat = status)
           if (status > 0) goto 9000
           
           call dsyev('V', & ! both ew's and ev's 
                'U', & ! upper triangle of A
                n, w%A, n, & ! data about A
                w%ew, workquery, -1, & 
                lapack_status)
        else
           allocate( w%iwork(5*n), stat = status )
           if (status > 0) goto 9000
           allocate( w%ifail(n), stat = status ) 
           if (status > 0) goto 9000
           
           ! make a workspace query to dsyevx
           call dsyevx( 'V',& ! get both ew's and ev's
                     'I',& ! just the numbered eigenvalues
                     'U',& ! upper triangle of A
                      n, w%A, n, & 
                      1.0, 1.0, & ! not used for RANGE = 'I'
                      1, 1, & ! only find the first eigenpair
                      0.5, & ! abstol for the eigensolver
                      eigsout, & ! total number of eigs found
                      1.0, 1.0, & ! the eigenvalue and eigenvector
                      n, & ! ldz (the eigenvector array)
                      workquery, -1, w%iwork, &  ! workspace
                      w%ifail, & ! array containing indicies of non-converging ews
                      lapack_status)
           if (lapack_status > 0) goto 9020
        end if
        lwork = int(workquery(1))
        deallocate(workquery)
        allocate( w%work(lwork), stat = status )
        if (status > 0) goto 9000

        return
        
9000    continue
        ! Allocation errors : min_eig_symm
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''min_eig_symm'': ',&
                'not enough memory.' 
        end if
        
        return

9020    continue
        ! Error return from lapack routine
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for lapack subroutine: ',&
                'not enough memory.' 
        end if
        return

      end subroutine setup_workspace_min_eig_symm
      
      subroutine setup_workspace_max_eig(n,m,w,control,status)
        integer, intent(in) :: n, m 
        type( max_eig_work) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        real(wp), allocatable :: workquery(:)
        integer :: lapack_status, lwork
        
        allocate( w%alphaR(2*n), stat = status)
        if (status > 0) goto 9000
        allocate( w%alphaI(2*n), stat = status)
        if (status > 0) goto 9000
        allocate( w%beta(2*n),   stat = status)
        if (status > 0) goto 9000
        allocate( w%vr(2*n,2*n), stat = status)
        if (status > 0) goto 9000
        allocate( w%ew_array(2*n), stat = status)
        if (status > 0) goto 9000
        allocate(workquery(1),stat = status)
        if (status > 0) goto 9000
        ! make a workspace query to dggev
        call dggev('N', & ! No left eigenvectors
             'V', &! Yes right eigenvectors
             2*n, 1.0, 2*n, 1.0, 2*n, &
             1.0, 0.1, 0.1, & ! eigenvalue data
             0.1, 2*n, & ! not referenced
             0.1, 2*n, & ! right eigenvectors
             workquery, -1, lapack_status)
        if (lapack_status > 0) goto 9020
        lwork = int(workquery(1))
        deallocate(workquery)
        allocate( w%work(lwork), stat = status)
        if (status > 0) goto 9000
        allocate( w%nullindex(2*n), stat = status)
        if (status > 0) goto 9000
        allocate( w%vecisreal(2*n), stat = status)
        if (status > 0) goto 9000

        return
        
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''AINT_tr'': ',&
                'not enough memory.' 
        end if
        
        return

9020    continue
        ! Error return from lapack routine
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for lapack subroutine: ',&
                'not enough memory.' 
        end if
        return

      end subroutine setup_workspace_max_eig


      subroutine setup_workspace_solve_general(n, m, w, control, status)
        integer, intent(in) :: n, m 
        type( solve_general_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        
        allocate( w%A(n,n), stat = status)
        if (status > 0) goto 9000
        allocate( w%ipiv(n), stat = status)
        if (status > 0) goto 9000
        
        return
        
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''solve_general'': ',&
                'not enough memory.' 
        end if
        
        return

      end subroutine setup_workspace_solve_general

      subroutine setup_workspace_solve_dtrs(n,m,w,control,status)
        integer, intent(in) :: n,m
        type( solve_dtrs_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        
        allocate(w%A(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%ev(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%v(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%v_trans(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%ew(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%d_trans(n),stat = status)
        if (status > 0) goto 9000

        call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,control,status)
        if (status > 0) goto 9010
        
        return
                
9000    continue
        ! Allocation errors : solve_dtrs
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''solve_dtrs'': ',&
                'not enough memory.' 
        end if
        
        return
        
9010    continue  
        ! errors : solve_dtrs
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Called from routine solve_dtrs' 
        end if
        
        return

      end subroutine setup_workspace_solve_dtrs
      
      subroutine setup_workspace_all_eig_symm(n,m,w,control,status)
        integer, intent(in) :: n,m
        type( all_eig_symm_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status

        real(wp), allocatable :: workquery(:)
        real(wp) :: A,ew
        integer :: lapack_status, lwork, eigsout

        lapack_status = 0
        A = 1.0_wp
        ew = 1.0_wp

        allocate(workquery(1))
        call dsyev('V', & ! both ew's and ev's 
             'U', & ! upper triangle of A
             n, A, n, & ! data about A
             ew, workquery, -1, & 
             lapack_status)  
        if (lapack_status .ne. 0) goto 9000

        lwork = int(workquery(1))
        deallocate(workquery)
        allocate( w%work(lwork), stat = status )
        if (status > 0) goto 8000
        
        return
        
8000    continue 
        ! Allocation errors : all_eig_sym
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''all_eig_symm'': ',&
                'not enough memory.' 
        end if

9000    continue
        ! lapack error
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''solve_dtrs'': ',&
                'not enough memory.' 
        end if
        
        return

      end subroutine setup_workspace_all_eig_symm

      
      subroutine setup_workspace_more_sorensen(n,m,w,control,status)
        integer, intent(in) :: n,m
        type( more_sorensen_work ) :: w
        type( NLLS_control_type ), intent(in) :: control
        integer, intent(out) :: status
        allocate(w%A(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%LtL(n,n),stat = status)
        if (status > 0) goto 9000
        allocate(w%v(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%q(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%y1(n),stat = status)
        if (status > 0) goto 9000
        allocate(w%AplusSigma(n,n),stat = status)
        if (status > 0) goto 9000

        call setup_workspace_min_eig_symm(n,m,w%min_eig_symm_ws,control,status)
        if (status > 0) goto 9010
        
        return
        
9000    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a,a)') &
                'Error allocating array for subroutine ''more_sorensen'': ',&
                'not enough memory.' 
        end if
        
        return
        
9010    continue
        ! Allocation errors : dogleg
        if (control%print_level >= 0) then
           write(control%error,'(a)') &
                'Called from subroutine ''dogleg'': '
        end if

        return


      end subroutine setup_workspace_more_sorensen
      

end module nlls_module

!!
!! Below are the modules that are borrowed from Galahad...
!!

! THIS VERSION: RAL_NLLS 1.0 - 22/12/2015 AT 14:15 GMT.

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*     SYMBOLS  M O D U L E    *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package SYMBOLS, December 22nd, 2015

   MODULE RAL_NLLS_SYMBOLS

!  This module provides the list of all symbolic names that are common to
!  the RAL_NLLS modules. It is just intended as a dictionary for use in other
!  modules.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      IMPLICIT NONE
      PUBLIC :: SYMBOLS_status

!  exit conditions (0 to -99; others will be package specific)

      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_ok                      = 0
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_allocate          = - 1
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_deallocate        = - 2
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_restrictions      = - 3
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_bad_bounds        = - 4
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_primal_infeasible = - 5
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_dual_infeasible   = - 6
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unbounded         = - 7
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_no_center         = - 8
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_analysis          = - 9
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_factorization     = - 10
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_solve             = - 11
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_analysis      = - 12
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_factorization = - 13
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_solve         = - 14
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_preconditioner    = - 15
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ill_conditioned   = - 16
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_tiny_step         = - 17
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_max_iterations    = - 18
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_time_limit        = - 19
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cpu_limit         =         &
                                    RAL_NLLS_error_time_limit
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_inertia           = - 20
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_file              = - 21
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_io                = - 22
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_upper_entry       = - 23
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_sort              = - 24
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_input_status      = - 25
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unknown_solver    = - 26
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_not_yet_implemented     = - 27
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qp_solve          = - 28
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_unavailable_option      = - 29
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_warning_on_boundary     = - 30
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_call_order        = - 31
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_integer_ws        = - 32
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_real_ws           = - 33
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_pardiso           = - 34
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_wsmp              = - 35
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc64              = - 36
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc77              = - 37
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_lapack            = - 38
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_permutation       = - 39
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_alter_diagonal    = - 40
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_access_pivots     = - 41
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_access_pert       = - 42
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_direct_access     = - 43
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_f_min             = - 44
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unknown_precond   = - 45
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_schur_complement  = - 46
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_technical         = - 50
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_reformat          = - 52
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ah_unordered      = - 53
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_y_unallocated     = - 54
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_z_unallocated     = - 55
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_scale             = - 61
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_presolve          = - 62
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpa               = - 63
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpb               = - 64
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpc               = - 65
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cqp               = - 66
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_dqp               = - 67
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc61              = - 69
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc68              = - 70
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_metis             = - 71
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_spral             = - 72
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_warning_repeated_entry  = - 73
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_rif               = - 74
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ls28              = - 75
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ls29              = - 76
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cutest            = - 77
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_evaluation        = - 78
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_optional          = - 79
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mi35              = - 80
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_spqr              = - 81
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_alive             = - 82
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ccqp              = - 83

   CONTAINS

!-*-  R A L _ N L L S -  S Y M B O L S _ S T A T U S   S U B R O U T I N E  -*-

     SUBROUTINE SYMBOLS_status( status, out, prefix, routine )

!  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!  Print details of return codes

!  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: status, out
     CHARACTER ( LEN = * ) :: prefix
     CHARACTER ( LEN = * ) :: routine

     SELECT CASE ( status )
     CASE( RAL_NLLS_ok )

     CASE( RAL_NLLS_error_allocate )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   allocation error' )" )                                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_deallocate )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   deallocation error' )" )                                     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_restrictions )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
     &        '   one or more restrictions violated' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_bad_bounds )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem bounds are infeasible' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_primal_infeasible )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem appears to be (locally) infeasible' )" )         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_dual_infeasible )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem appears to be (locally) dual infeasible' )" )    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unbounded )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
     &        '   the problem appears to be unbounded from below' )" )         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_no_center )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the feasible region has no analytic center' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_analysis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the analysis phase failed' )" )                              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_factorization )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the factorization failed' )" )                               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the solve phase failed' )" )                                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_analysis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &        '   the unsymmetric analysis phase failed' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_factorization )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the unsymmetric factorization failed' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the unsymmetric solve phase failed' )" )                     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_preconditioner )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the provided preconditioner is flawed' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ill_conditioned )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A, '   the',   &
      &       ' problem is too ill conditioned to make further progress' )" )  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_tiny_step )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A, '   the',   &
      &       ' computed step is too small to make further progress' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_max_iterations )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the iteration limit has been exceeded' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_time_limit )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the time limit has been exceeded' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_inertia )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the preconditioner has the wrong inertia' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_file )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a file-handling error occurred' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_io )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an input/output error occurred' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_upper_entry )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   there is a matrix entry in the upper triangle' )" )          &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_sort )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occurred when sorting' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_input_status )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   error with input status' )" )                                &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unknown_solver )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested solver is not known' )" )                      &
         prefix, routine, prefix
     CASE ( RAL_NLLS_not_yet_implemented )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested option has not yet been implemented' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qp_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   QP solver failed: check the QP solver status' )" )           &
         prefix, routine, prefix
     CASE( RAL_NLLS_unavailable_option )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested option is unavailable' )" )                    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_call_order )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the routine has been called out of order' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_integer_ws )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   insufficient integer workspace for factorization' )" )       &
         prefix, routine, prefix, status
     CASE( RAL_NLLS_error_real_ws )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   insufficient real workspace for factorization' )" )          &
         prefix, routine, prefix, status
     CASE( RAL_NLLS_error_pardiso )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   PARDISO failure: check its return status' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_wsmp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   WSMP failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc64 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC64 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc77 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC77 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_lapack )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LAPACK failure: check its return status' )" )                &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_permutation )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   input order is not a permutation' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_alter_diagonal )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to alter diagonals with this solver' )" )   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_access_pivots )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to access the pivots with this solver' )" ) &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_access_pert )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to access the perturbations made with',     &
      &       ' this solver' )" ) prefix, routine, prefix
     CASE( RAL_NLLS_error_direct_access )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a direct-access file error occurred' )" )                    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_f_min )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the objective value is too small' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unknown_precond )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested preconditioner is not known' )" )              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_schur_complement )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested Schur complement is too large' )" )            &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_technical )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a technical error has occurred within the linear solver' )" )&
         prefix, routine, prefix
     CASE( RAL_NLLS_error_reformat )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the matrix storage format has been changed' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ah_unordered )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the matrix storage format should have been changed' )" )     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_y_unallocated )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   Y, Y_l or Y_u has not been allocated' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_z_unallocated )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   Z, Z_l or Z_u has not been allocated' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_scale )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when scaling the problem' )" )              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_presolve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when pre/postsolving the problem' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpa )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPA' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpb )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPB' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpc )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPC' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_cqp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling CQP' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_dqp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling DQP' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc61 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC61 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc68 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC68 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_metis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   METIS failure: check its return status' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_spral )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   SPRAL failure: check its return status' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ls28 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LS28 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ls29 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LS29 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mi35 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MI35 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_rif )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   RIF failure: check its return status' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_warning_repeated_entry )
       WRITE( out, "( /, A,  ' Warning return from ', A, ' -', /, A,           &
      &       '   the input matrix contains repeated entries' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_cutest )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   CUTEst evaluation error: check its return status' )" )       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_evaluation )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   function evaluation error: check its return status' )" )     &
         prefix, routine, prefix
     CASE DEFAULT
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   status = ', I0 )" ) prefix, routine, prefix, status
     END SELECT

     RETURN

!  End of subroutine SYMBOLS_status

     END SUBROUTINE SYMBOLS_status

   END MODULE RAL_NLLS_SYMBOLS


! THIS VERSION: RAL_NLLS 1.0 - 22/12/2015 AT 14:15 GMT.

!-*-*-*-*-*-*-*-  R A L _ N L L S _ D T R S  double  M O D U L E  *-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package TRS, December 22nd, 2015

   MODULE RAL_NLLS_DTRS_double

!       -----------------------------------------------
!      |                                               |
!      | Solve the trust-region subproblem             |
!      |                                               |
!      |    minimize     1/2 <x, H x> + <c, x> + f     |
!      |    subject to      ||x||_2 <= radius          |
!      |    or              ||x||_2  = radius          |
!      |                                               |
!      | where H is diagonal                           |
!      |                                               |
!       -----------------------------------------------

      USE RAL_NLLS_SYMBOLS

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: DTRS_initialize, DTRS_solve

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: history_max = 100
      INTEGER, PARAMETER :: max_degree = 3
      INTEGER, PARAMETER :: out = 6
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: point4 = 0.4_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: four = 4.0_wp
      REAL ( KIND = wp ), PARAMETER :: six = 6.0_wp
      REAL ( KIND = wp ), PARAMETER :: sixth = one / six
      REAL ( KIND = wp ), PARAMETER :: onethird = one / three
      REAL ( KIND = wp ), PARAMETER :: twothirds = two /three
      REAL ( KIND = wp ), PARAMETER :: onesixth = one / six
      REAL ( KIND = wp ), PARAMETER :: threequarters = 0.75_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: twentyfour = 24.0_wp
      REAL ( KIND = wp ), PARAMETER :: infinity = half * HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: teneps = ten * epsmch
      REAL ( KIND = wp ), PARAMETER :: pi = 3.1415926535897931_wp
      REAL ( KIND = wp ), PARAMETER :: magic = 2.0943951023931953_wp  !! 2 pi/3
      REAL ( KIND = wp ), PARAMETER :: roots_tol = teneps
      LOGICAL :: roots_debug = .FALSE.

!--------------------------
!  Derived type definitions
!--------------------------

!  - - - - - - - - - - - - - - - - - - - - - - -
!   control derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_control_type

!  unit for error messages

        INTEGER :: error = 6

!  unit for monitor output

        INTEGER :: out = 6

!  unit to write problem data into file problem_file

        INTEGER :: problem = 0

!  controls level of diagnostic output

        INTEGER :: print_level = 0

!  maximum degree of Taylor approximant allowed

        INTEGER :: taylor_max_degree = 3

!  lower and upper bounds on the multiplier, if known

        REAL ( KIND = wp ) :: lower = - half * HUGE( one )
        REAL ( KIND = wp ) :: upper =  HUGE( one )

!  stop when | ||x|| - radius | <=
!     max( stop_normal * radius, stop_absolute_normal )

        REAL ( KIND = wp ) :: stop_normal = epsmch
        REAL ( KIND = wp ) :: stop_absolute_normal = epsmch

!  is the solution is REQUIRED to lie on the boundary (i.e., is the constraint
!  an equality)?

        LOGICAL :: equality_problem= .FALSE.

!  name of file into which to write problem data

        CHARACTER ( LEN = 30 ) :: problem_file =                               &
         'trs_problem.data' // REPEAT( ' ', 14 )

!  all output lines will be prefixed by
!    prefix(2:LEN(TRIM(%prefix))-1)
!  where prefix contains the required string enclosed in quotes,
!  e.g. "string" or 'string'

        CHARACTER ( LEN = 30 ) :: prefix  = '""                            '

      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - - -
!   history derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_history_type

!  value of lambda

        REAL ( KIND = wp ) :: lambda = zero

!  corresponding value of ||x(lambda)||_M

        REAL ( KIND = wp ) :: x_norm = zero
      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - -
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_inform_type

!   reported return status:
!      0 the solution has been found
!     -3 n and/or Delta is not positive
!    -16 ill-conditioning has prevented furthr progress

        INTEGER :: status = 0

!  the number of (||x||_M,lambda) pairs in the history

        INTEGER :: len_history = 0

!  the value of the quadratic function

        REAL ( KIND = wp ) :: obj = HUGE( one )

!  the M-norm of x, ||x||_M

        REAL ( KIND = wp ) :: x_norm = zero

!  the Lagrange multiplier corresponding to the trust-region constraint

        REAL ( KIND = wp ) :: multiplier = zero

!  a lower bound max(0,-lambda_1), where lambda_1 is the left-most
!  eigenvalue of (H,M)

        REAL ( KIND = wp ) :: pole = zero

!  has the hard case occurred?

        LOGICAL :: hard_case = .FALSE.

!  history information

        TYPE ( DTRS_history_type ), DIMENSION( history_max ) :: history
      END TYPE

!  interface to LAPACK: eigenvalues of a Hessenberg matrix

      INTERFACE HSEQR
        SUBROUTINE SHSEQR( job, compz, n, ilo, ihi, H, ldh,  WR, WI, Z, ldz,   &
                           WORK, lwork, info )
        INTEGER, INTENT( IN ) :: ihi, ilo, ldh, ldz, lwork, n
        INTEGER, INTENT( OUT ) :: info
        CHARACTER ( LEN = 1 ), INTENT( IN ) :: compz, job
        REAL, INTENT( INOUT ) :: H( ldh, * ), Z( ldz, * )
        REAL, INTENT( OUT ) :: WI( * ), WR( * ), WORK( * )
        END SUBROUTINE SHSEQR

        SUBROUTINE DHSEQR( job, compz, n, ilo, ihi, H, ldh,  WR, WI, Z, ldz,   &
                           WORK, lwork, info )
        INTEGER, INTENT( IN ) :: ihi, ilo, ldh, ldz, lwork, n
        INTEGER, INTENT( OUT ) :: info
        CHARACTER ( LEN = 1 ), INTENT( IN ) :: compz, job
        DOUBLE PRECISION, INTENT( INOUT ) :: H( ldh, * ), Z( ldz, * )
        DOUBLE PRECISION, INTENT( OUT ) :: WI( * ), WR( * ), WORK( * )
        END SUBROUTINE DHSEQR
      END INTERFACE HSEQR

!  interface to BLAS: two norm

     INTERFACE NRM2

       FUNCTION SNRM2( n, X, incx )
       REAL :: SNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION SNRM2

       FUNCTION DNRM2( n, X, incx )
       DOUBLE PRECISION :: DNRM2
       INTEGER, INTENT( IN ) :: n, incx
       DOUBLE PRECISION, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION DNRM2

     END INTERFACE

    CONTAINS

!-*-*-*-*-*-*-  D T R S _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE DTRS_initialize( control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  .  Set initial values for the TRS control parameters  .
!
!  Arguments:
!  =========
!
!   control  a structure containing control information. See DTRS_control_type
!   data     private internal data
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------

      TYPE ( DTRS_CONTROL_TYPE ), INTENT( OUT ) :: control
      TYPE ( DTRS_inform_type ), INTENT( OUT ) :: inform

      inform%status = RAL_NLLS_ok

!  Set initial control parameter values

      control%stop_normal = epsmch ** 0.75
      control%stop_absolute_normal = epsmch ** 0.75

      RETURN

!  End of subroutine DTRS_initialize

      END SUBROUTINE DTRS_initialize

!-*-*-*-*-*-*-*-*-  D T R S _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*-

      SUBROUTINE DTRS_solve( n, radius, f, C, H, X, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the trust-region subproblem
!
!      minimize     1/2 <x, H x> + <c, x> + f
!      subject to    ||x||_2 <= radius  or ||x||_2 = radius
!
!  where H is diagonal, using a secular iteration
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   radius - the trust-region radius
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a vector of values for the diagonal matrix H
!
!   X - the required solution vector x
!
!   control - a structure containing control information. See DTRS_control_type
!
!   inform - a structure containing information. See DTRS_inform_type
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: radius
      REAL ( KIND = wp ), INTENT( IN ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C, H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( DTRS_control_type ), INTENT( IN ) :: control
      TYPE ( DTRS_inform_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, it, out, nroots, print_level
      INTEGER :: max_order, n_lambda, i_hard
      REAL ( KIND = wp ) :: lambda, lambda_l, lambda_u, delta_lambda
      REAL ( KIND = wp ) :: alpha, utx, distx
      REAL ( KIND = wp ) :: c_norm, c_norm_over_radius, v_norm2, w_norm2
      REAL ( KIND = wp ) :: beta, z_norm2, root1, root2, root3
      REAL ( KIND = wp ) :: lambda_min, lambda_max, lambda_plus, c2
      REAL ( KIND = wp ) :: a_0, a_1, a_2, a_3, a_max
      REAL ( KIND = wp ), DIMENSION( 3 ) :: lambda_new
      REAL ( KIND = wp ), DIMENSION( 0 : max_degree ) :: x_norm2, pi_beta
      LOGICAL :: printi, printt, printd, problem_file_exists
      CHARACTER ( LEN = 1 ) :: region

!  prefix for all output

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      IF ( LEN( TRIM( control%prefix ) ) > 2 )                                 &
        prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  output problem data
      IF ( control%problem > 0 ) THEN
        INQUIRE( FILE = control%problem_file, EXIST = problem_file_exists )
        IF ( problem_file_exists ) THEN
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'OLD' )
          REWIND control%problem
        ELSE
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'NEW' )
        END IF
        WRITE( control%problem, * ) n, COUNT( C( : n ) /= zero ),              &
          COUNT( H( : n ) /= zero )
        WRITE( control%problem, * ) radius, f
        DO i = 1, n
          IF ( C( i ) /= zero ) WRITE( control%problem, * ) i, C( i )
        END DO
        DO i = 1, n
          IF ( H( i ) /= zero ) WRITE( control%problem, * ) i, i, H( i )
        END DO
        CLOSE( control%problem )
      END IF

!  set initial values

      X = zero ; inform%x_norm = zero ; inform%obj = f
      inform%hard_case = .FALSE.
      delta_lambda = zero

!  record desired output level

      out = control%out
      print_level = control%print_level
      printi = out > 0 .AND. print_level > 0
      printt = out > 0 .AND. print_level > 1
      printd = out > 0 .AND. print_level > 2

!  compute the two-norm of c and the extreme eigenvalues of H

      c_norm = TWO_NORM( C )
      lambda_min = MINVAL( H( : n ) )
      lambda_max = MAXVAL( H( : n ) )

      IF ( printt ) WRITE( out, "( A, ' ||c|| = ', ES10.4, ', ||H|| = ',       &
     &                             ES10.4, ', lambda_min = ', ES11.4 )" )      &
          prefix, c_norm, MAXVAL( ABS( H( : n ) ) ), lambda_min

      region = 'L'
      IF ( printt )                                                            &
        WRITE( out, "( A, 4X, 28( '-' ), ' phase two ', 28( '-' ) )" ) prefix
      IF ( printi ) WRITE( out, 2030 ) prefix

!  check for the trivial case

      IF ( c_norm == zero .AND. lambda_min >= zero ) THEN
        IF (  control%equality_problem ) THEN
          DO i = 1, n
            IF ( H( i ) == lambda_min ) THEN
              i_hard = i
              EXIT
            END IF
          END DO
          X( i_hard ) = one / radius
          inform%x_norm = radius
          inform%obj = f + lambda_min * radius ** 2
          lambda = - lambda_min
        ELSE
          lambda = zero
        END IF
        IF ( printi ) THEN
          WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,             &
          0, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
          WRITE( out, "( A,                                                    &
       &    ' Normal stopping criteria satisfied' )" ) prefix
        END IF
        inform%status = RAL_NLLS_ok
        GO TO 900
      END IF

!  construct values lambda_l and lambda_u for which lambda_l <= lambda_optimal
!   <= lambda_u, and ensure that all iterates satisfy lambda_l <= lambda
!   <= lambda_u

      c_norm_over_radius = c_norm / radius
      IF ( control%equality_problem ) THEN
        lambda_l = MAX( control%lower,  - lambda_min,                          &
                        c_norm_over_radius - lambda_max )
        lambda_u = MIN( control%upper,                                         &
                        c_norm_over_radius - lambda_min )
      ELSE
        lambda_l = MAX( control%lower,  zero, - lambda_min,                    &
                        c_norm_over_radius - lambda_max )
        lambda_u = MIN( control%upper,                                         &
                        MAX( zero, c_norm_over_radius - lambda_min ) )
      END IF
      lambda = lambda_l

!  check for the "hard case"

      IF ( lambda == - lambda_min ) THEN
        c2 = zero
        inform%hard_case = .TRUE.
        DO i = 1, n
          IF ( H( i ) == lambda_min ) THEN
            IF ( ABS( C( i ) ) > epsmch * c_norm ) THEN
              inform%hard_case = .FALSE.
              c2 = c2 + C( i ) ** 2
            ELSE
              i_hard = i
            END IF
          END IF
        END DO

!  the hard case may occur

        IF ( inform%hard_case ) THEN
          DO i = 1, n
            IF ( H( i ) /= lambda_min ) THEN
              X( i )  = - C( i ) / ( H( i ) + lambda )
            ELSE
              X( i ) = zero
            END IF
          END DO
          inform%x_norm = TWO_NORM( X )

!  the hard case does occur

          IF ( inform%x_norm <= radius ) THEN
            IF ( inform%x_norm < radius ) THEN

!  compute the step alpha so that X + alpha E_i_hard lies on the trust-region
!  boundary and gives the smaller value of q

              utx = X( i_hard ) / radius
              distx = ( radius - inform%x_norm ) *                             &
                ( ( radius + inform%x_norm ) / radius )
              alpha = sign( distx / ( abs( utx ) +                             &
                            sqrt( utx ** 2 + distx / radius ) ), utx )

!  record the optimal values

              X( i_hard ) = X( i_hard ) + alpha
            END IF
            inform%x_norm = TWO_NORM( X )
            inform%obj =                                                       &
                f + half * ( DOT_PRODUCT( C, X ) - lambda * radius ** 2 )
            IF ( printi ) THEN
              WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,         &
              0, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
              WRITE( out, "( A,                                                &
           &    ' Hard-case stopping criteria satisfied' )" ) prefix
            END IF
            inform%status = RAL_NLLS_ok
            GO TO 900

!  the hard case didn't occur after all

          ELSE
            inform%hard_case = .FALSE.

!  compute the first derivative of ||x|(lambda)||^2 - radius^2

            w_norm2 = zero
            DO i = 1, n
              IF ( H( i ) /= lambda_min )                                      &
                w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
            END DO
            x_norm2( 1 ) = - two * w_norm2

!  compute the Newton correction

            lambda =                                                           &
              lambda + ( inform%x_norm ** 2 - radius ** 2 ) / x_norm2( 1 )
            lambda_l = MAX( lambda_l, lambda )
          END IF

!  there is a singularity at lambda. Compute the point for which the
!  sum of squares of the singular terms is equal to radius^2

        ELSE
          lambda = lambda + SQRT( c2 ) / radius
          lambda_l = MAX( lambda_l, lambda )
        END IF
      END IF

!  the iterates will all be in the L region. Prepare for the main loop

      it = 0
      max_order = MAX( 1, MIN( max_degree, control%taylor_max_degree ) )

!  start the main loop

      DO
        it = it + 1

!  if H(lambda) is positive definite, solve  H(lambda) x = - c

        DO i = 1, n
          X( i )  = - C( i ) / ( H( i ) + lambda )
        END DO

!  compute the two-norm of x

        inform%x_norm = TWO_NORM( X )
        x_norm2( 0 ) = inform%x_norm ** 2

!  if the Newton step lies within the trust region, exit

        IF ( lambda == zero .AND. inform%x_norm <= radius ) THEN
          inform%obj = f + half * DOT_PRODUCT( C, X )
          inform%status = RAL_NLLS_ok
          region = 'L'
          IF ( printi ) THEN
            WRITE( out, "( A, A2, I4, 2ES22.15 )" ) prefix,                    &
              region, it, inform%x_norm - radius, lambda
            WRITE( out, "( A, ' Interior stopping criteria satisfied')" ) prefix
          END IF
          GO TO 900
        END IF

!  the current estimate gives a good approximation to the required root

        IF ( ABS( inform%x_norm - radius ) <=                                  &
               MAX( control%stop_normal * radius,                              &
                    control%stop_absolute_normal ) ) THEN
          IF ( inform%x_norm > radius ) THEN
            lambda_l = MAX( lambda_l, lambda )
          ELSE
            region = 'G'
            lambda_u = MIN( lambda_u, lambda )
          END IF
          IF ( printt .AND. it > 1 ) WRITE( out, 2030 ) prefix
          IF ( printi ) THEN
            WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,           &
            it, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
            WRITE( out, "( A,                                                  &
         &    ' Normal stopping criteria satisfied' )" ) prefix
          END IF
          inform%status = RAL_NLLS_ok
          EXIT
        END IF

        lambda_l = MAX( lambda_l, lambda )

!  debug printing

        IF ( printd ) THEN
          WRITE( out, "( A, 8X, 'lambda', 13X, 'x_norm', 15X, 'radius' )" )    &
            prefix
          WRITE( out, "( A, 3ES20.12 )") prefix, lambda, inform%x_norm, radius
        ELSE IF ( printi ) THEN
          IF ( printt .AND. it > 1 ) WRITE( out, 2030 ) prefix
          WRITE( out, "( A, A2, I4, 3ES22.15 )" ) prefix, region, it,          &
            ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
        END IF

!  record, for the future, values of lambda which give small ||x||

        IF ( inform%len_history < history_max ) THEN
          inform%len_history = inform%len_history + 1
          inform%history( inform%len_history )%lambda = lambda
          inform%history( inform%len_history )%x_norm = inform%x_norm
        END IF

!  a lambda in L has been found. It is now simply a matter of applying
!  a variety of Taylor-series-based methods starting from this lambda

!  precaution against rounding producing lambda outside L

        IF ( lambda > lambda_u ) THEN
          inform%status = RAL_NLLS_error_ill_conditioned
          EXIT
        END IF

!  compute first derivatives of x^T M x

!  form ||w||^2 = x^T H^-1(lambda) x

        w_norm2 = zero
        DO i = 1, n
          w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
        END DO

!  compute the first derivative of x_norm2 = x^T M x

        x_norm2( 1 ) = - two * w_norm2

!  compute pi_beta = ||x||^beta and its first derivative when beta = - 1

        beta = - one
        CALL DTRS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )

!  compute the Newton correction (for beta = - 1)

        delta_lambda = - ( pi_beta( 0 ) - ( radius ) ** beta ) / pi_beta( 1 )

        n_lambda = 1
        lambda_new( n_lambda ) = lambda + delta_lambda

        IF ( max_order >= 3 ) THEN

!  compute the second derivative of x^T x

          z_norm2 = zero
          DO i = 1, n
            z_norm2 = z_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 4
          END DO
          x_norm2( 2 ) = six * z_norm2

!  compute the third derivatives of x^T x

          v_norm2 = zero
          DO i = 1, n
            v_norm2 = v_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 5
          END DO
          x_norm2( 3 ) = - twentyfour * v_norm2

!  compute pi_beta = ||x||^beta and its derivatives when beta = 2

          beta = two
          CALL DTRS_pi_derivs( max_order, beta, x_norm2( : max_order ),        &
                               pi_beta( : max_order ) )

!  compute the "cubic Taylor approximaton" step (beta = 2)

          a_0 = pi_beta( 0 ) - ( radius ) ** beta
          a_1 = pi_beta( 1 )
          a_2 = half * pi_beta( 2 )
          a_3 = sixth * pi_beta( 3 )
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max
            a_2 = a_2 / a_max ; a_3 = a_3 / a_max
          END IF
          CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,             &
                            root1, root2, root3, roots_debug )
          n_lambda = n_lambda + 1
          IF ( nroots == 3 ) THEN
            lambda_new( n_lambda ) = lambda + root3
          ELSE
            lambda_new( n_lambda ) = lambda + root1
          END IF

!  compute pi_beta = ||x||^beta and its derivatives when beta = - 0.4

          beta = - point4
          CALL DTRS_pi_derivs( max_order, beta, x_norm2( : max_order ),        &
                               pi_beta( : max_order ) )

!  compute the "cubic Taylor approximaton" step (beta = - 0.4)

          a_0 = pi_beta( 0 ) - ( radius ) ** beta
          a_1 = pi_beta( 1 )
          a_2 = half * pi_beta( 2 )
          a_3 = sixth * pi_beta( 3 )
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max
            a_2 = a_2 / a_max ; a_3 = a_3 / a_max
          END IF
          CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,             &
                            root1, root2, root3, roots_debug )
          n_lambda = n_lambda + 1
          IF ( nroots == 3 ) THEN
            lambda_new( n_lambda ) = lambda + root3
          ELSE
            lambda_new( n_lambda ) = lambda + root1
          END IF
        END IF

!  record all of the estimates of the optimal lambda

        IF ( printd ) THEN
          WRITE( out, "( A, ' lambda_t (', I1, ')', 3ES20.13 )" )              &
            prefix, MAXLOC( lambda_new( : n_lambda ) ),                        &
            lambda_new( : MIN( 3, n_lambda ) )
          IF ( n_lambda > 3 ) WRITE( out, "( A, 13X, 3ES20.13 )" )             &
            prefix, lambda_new( 4 : MIN( 6, n_lambda ) )
        END IF

!  compute the best Taylor improvement

        lambda_plus = MAXVAL( lambda_new( : n_lambda ) )
        delta_lambda = lambda_plus - lambda
        lambda = lambda_plus

!  improve the lower bound if possible

        lambda_l = MAX( lambda_l, lambda_plus )

!  check that the best Taylor improvement is significant

        IF ( ABS( delta_lambda ) < epsmch * MAX( one, ABS( lambda ) ) ) THEN
          IF ( printi ) WRITE( out, "( A, ' normal exit with no ',             &
         &                     'significant Taylor improvement' )" ) prefix
          inform%status = RAL_NLLS_ok
          EXIT
        END IF

!  End of main iteration loop

      END DO

!  Record the optimal obective value

      inform%obj = f + half * ( DOT_PRODUCT( C, X ) - lambda * x_norm2( 0 ) )

!  ----
!  Exit
!  ----

 900  CONTINUE
      inform%multiplier = lambda
      inform%pole = MAX( zero, - lambda_min )
      RETURN

! Non-executable statements

 2030 FORMAT( A, '    it     ||x||-radius             lambda ',                &
                 '              d_lambda' )

!  End of subroutine DTRS_solve

      END SUBROUTINE DTRS_solve

!-*-*-*-*-*-*-  D T R S _ P I _ D E R I V S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE DTRS_pi_derivs( max_order, beta, x_norm2, pi_beta )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute pi_beta = ||x||^beta and its derivatives
!
!  Arguments:
!  =========
!
!  Input -
!   max_order - maximum order of derivative
!   beta - power
!   x_norm2 - (0) value of ||x||^2,
!             (i) ith derivative of ||x||^2, i = 1, max_order
!  Output -
!   pi_beta - (0) value of ||x||^beta,
!             (i) ith derivative of ||x||^beta, i = 1, max_order
!
!  Extracted wholesale from module RAL_NLLS_RQS
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: max_order
      REAL ( KIND = wp ), INTENT( IN ) :: beta, x_norm2( 0 : max_order )
      REAL ( KIND = wp ), INTENT( OUT ) :: pi_beta( 0 : max_order )

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      REAL ( KIND = wp ) :: hbeta

      hbeta = half * beta
      pi_beta( 0 ) = x_norm2( 0 ) ** hbeta
      pi_beta( 1 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - one ) ) * x_norm2( 1 )
      IF ( max_order == 1 ) RETURN
      pi_beta( 2 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - two ) ) *             &
        ( ( hbeta - one ) * x_norm2( 1 ) ** 2 + x_norm2( 0 ) * x_norm2( 2 ) )
      IF ( max_order == 2 ) RETURN
      pi_beta( 3 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - three ) ) *           &
        ( x_norm2( 3 ) * x_norm2( 0 ) ** 2 + ( hbeta - one ) *                 &
          ( three * x_norm2( 0 ) * x_norm2( 1 ) * x_norm2( 2 ) +               &
            ( hbeta - two ) * x_norm2( 1 ) ** 3 ) )

      RETURN

!  End of subroutine DTRS_pi_derivs

      END SUBROUTINE DTRS_pi_derivs

!-*-*-*-*-*  D T R S _ T H E T A _ D E R I V S   S U B R O U T I N E   *-*-*-*-

      SUBROUTINE DTRS_theta_derivs( max_order, beta, lambda, sigma,            &
                                     theta_beta )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute theta_beta = (lambda/sigma)^beta and its derivatives
!
!  Arguments:
!  =========
!
!  Input -
!   max_order - maximum order of derivative
!   beta - power
!   lambda, sigma - lambda and sigma
!  Output -
!   theta_beta - (0) value of (lambda/sigma)^beta,
!             (i) ith derivative of (lambda/sigma)^beta, i = 1, max_order
!
!  Extracted wholesale from module RAL_NLLS_RQS
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: max_order
      REAL ( KIND = wp ), INTENT( IN ) :: beta, lambda, sigma
      REAL ( KIND = wp ), INTENT( OUT ) :: theta_beta( 0 : max_order )

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      REAL ( KIND = wp ) :: los, oos

      los = lambda / sigma
      oos = one / sigma

!if ( los <= zero ) write( 6,*) ' l, s ',  lambda, sigma, beta
      theta_beta( 0 ) = los ** beta
      theta_beta( 1 ) = beta * ( los ** ( beta - one ) ) * oos
      IF ( max_order == 1 ) RETURN
      theta_beta( 2 ) = beta * ( los ** ( beta - two ) ) *                    &
                        ( beta - one ) * oos ** 2
      IF ( max_order == 2 ) RETURN
      theta_beta( 3 ) = beta * ( los ** ( beta - three ) ) *                  &
                        ( beta - one ) * ( beta - two ) * oos ** 3

      RETURN

!  End of subroutine DTRS_theta_derivs

      END SUBROUTINE DTRS_theta_derivs

!-*-*-*-*-*-   R O O T S _ q u a d r a t i c  S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2, debug )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the quadratic equation
!
!                   a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1 and a2 are real
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a2, a1, a0, tol
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2
      LOGICAL, INTENT( IN ) :: debug

!  Local variables

      REAL ( KIND = wp ) :: rhs, d, p, pprime

      rhs = tol * a1 * a1
      IF ( ABS( a0 * a2 ) > rhs ) THEN  !  really is quadratic
        root2 = a1 * a1 - four * a2 * a0
        IF ( ABS( root2 ) <= ( epsmch * a1 ) ** 2 ) THEN ! numerical double root
          nroots = 2 ; root1 = -  half * a1 / a2 ; root2 = root1
        ELSE IF ( root2 < zero ) THEN    ! complex not real roots
          nroots = 0 ; root1 = zero ; root2 = zero
        ELSE                             ! distint real roots
          d = - half * ( a1 + SIGN( SQRT( root2 ), a1 ) )
          nroots = 2 ; root1 = d / a2 ; root2 = a0 / d
          IF ( root1 > root2 ) THEN
            d = root1 ; root1 = root2 ; root2 = d
          END IF
        END IF
      ELSE IF ( a2 == zero ) THEN
        IF ( a1 == zero ) THEN
          IF ( a0 == zero ) THEN         ! the function is zero
            nroots = 1 ; root1 = zero ; root2 = zero
          ELSE                           ! the function is constant
            nroots = 0 ; root1 = zero ; root2 = zero
          END IF
        ELSE                             ! the function is linear
          nroots = 1 ; root1 = - a0 / a1 ; root2 = zero
        END IF
      ELSE                               ! very ill-conditioned quadratic
        nroots = 2
        IF ( - a1 / a2 > zero ) THEN
          root1 = zero ; root2 = - a1 / a2
        ELSE
          root1 = - a1 / a2 ; root2 = zero
        END IF
      END IF

!  perfom a Newton iteration to ensure that the roots are accurate

      IF ( nroots >= 1 ) THEN
        p = ( a2 * root1 + a1 ) * root1 + a0
        pprime = two * a2 * root1 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime
          root1 = root1 - p / pprime
          p = ( a2 * root1 + a1 ) * root1 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 1, root1, p
        IF ( nroots == 2 ) THEN
          p = ( a2 * root2 + a1 ) * root2 + a0
          pprime = two * a2 * root2 + a1
          IF ( pprime /= zero ) THEN
            IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime
            root2 = root2 - p / pprime
            p = ( a2 * root2 + a1 ) * root2 + a0
          END IF
          IF ( debug ) WRITE( out, 2010 ) 2, root2, p
        END IF
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quadratic = ', ES12.4,     &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quadratic = ', ES12.4 )


!  End of subroutine ROOTS_quadratic

      END SUBROUTINE ROOTS_quadratic

!-*-*-*-*-*-*-*-   R O O T S _ c u b i c  S U B R O U T I N E   -*-*-*-*-*-*-*-

      SUBROUTINE ROOTS_cubic( a0, a1, a2, a3, tol, nroots, root1, root2,       &
                              root3, debug )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the cubicc equation
!
!                a3 * x**3 + a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1, a2 and a3 are real
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a3, a2, a1, a0, tol
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2, root3
      LOGICAL, INTENT( IN ) :: debug

!  Local variables

      INTEGER :: info, nroots_q
      REAL ( KIND = wp ) :: a, b, c, d, e, f, p, q, s, t, w, x, y, z
      REAL ( KIND = wp ) :: c0, c1, c2, b0, b1, pprime, u1, u2
      REAL ( KIND = wp ) :: H( 3, 3 ), ER( 3 ), EI( 3 ), ZZ( 1, 3 ), WORK( 33 )

!  define method used:
!    1 = Nonweiler, 2 = Littlewood, 3 = Viete, other = companion matrix

      INTEGER, PARAMETER :: method = 1

!  Check to see if the quartic is actually a cubic

      IF ( a3 == zero ) THEN
        CALL ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2, debug )
        root3 = infinity
        RETURN
      END IF

!  Deflate the polnomial if the trailing coefficient is zero

      IF ( a0 == zero ) THEN
        root1 = zero
        CALL ROOTS_quadratic( a1, a2, a3, tol, nroots, root2, root3, debug )
        nroots = nroots + 1
        RETURN
      END IF

!  1. Use Nonweiler's method (CACM 11:4, 1968, pp269)

      IF ( method == 1 ) THEN
        c0 = a0 / a3
        c1 = a1 / a3
        c2 = a2 / a3

        s = c2 / three
        t = s * c2
        b = 0.5_wp * ( s * ( twothirds * t - c1 ) + c0 )
        t = ( t - c1 ) / three
        c = t * t * t ; d = b * b - c

! 1 real + 2 equal real or 2 complex roots

        IF ( d >= zero ) THEN
          d = ( SQRT( d ) + ABS( b ) ) ** onethird
          IF ( d /= zero ) then
            IF ( b > zero ) then
              b = - d
            ELSE
              b = d
            END IF
            c = t / b
          END IF
          d = SQRT( threequarters ) * ( b - c )
          b = b + c ; c = - 0.5 * b - s
          root1 = b - s
          IF ( d == zero ) THEN
            nroots = 3 ; root2 = c ; root3 = c
          ELSE
            nroots = 1
          END IF

! 3 real roots

        ELSE
          IF ( b == zero ) THEN
            d = twothirds * ATAN( one )
          ELSE
            d = ATAN( SQRT( - d ) / ABS( b ) ) / three
          END IF
          IF ( b < zero ) THEN
            b = two * SQRT( t )
          ELSE
            b = - two * SQRT( t )
          END IF
          c = COS( d ) * b
          t = - SQRT( threequarters ) * SIN( d ) * b - half * c
          d = - t - c - s ; c = c - s ; t = t - s
          IF ( ABS( c ) > ABS( t ) ) then
            root3 = c
          ELSE
            root3 = t
            t = c
          END IF
          IF ( ABS( d ) > ABS( t ) ) THEN
            root2 = d
          ELSE
            root2 = t
            t = d
          END IF
          root1 = t ; nroots = 3
        END IF

!  2. Use Littlewood's method

      ELSE IF ( method == 2 ) THEN
        c2 = a2 / ( three * a3 ) ; c1 = a1 / ( three * a3 ) ; c0 = a0 / a3
        x = c1 - c2 * c2
        y = c0 - c2* ( x + x + c1 )
        z = y ** 2 + four * x ** 3

!  there are three real roots

        IF ( z < zero ) THEN
          a = - two * SQRT( - x )
          b = y / ( a * x )
          y = ATAN2( SQRT( one - b ), SQRT( one + b ) ) * twothirds
          IF ( c2 < zero ) y = y + magic

!  calculate root which does not involve cancellation

          nroots = 1 ; root1 = a * COS( y ) - c2

!  there may be only one real root

        ELSE
          a = SQRT( z ) ; b = half * ( ABS( y ) + a ) ; c = b ** onethird
          IF ( c <= zero ) THEN
            nroots = 3 ; root1 = - c2 ; root2 = - c2 ; root3 = - c2
            GO TO 900
          ELSE
            nroots = 1
            c = c - ( c ** 3 - b ) / ( three * c * c )
            e = c * c + ABS( x )
            f = one / ( ( x / c ) ** 2 + e )
            IF ( x >= zero ) THEN
              x = e / c ; z = y * f
            ELSE
              x = a * f ; z = SIGN( one, y ) * e / c
            END IF
            IF ( z * c2 >= zero ) THEN
              root1 = - z - c2
            ELSE
              root2 = half * z - c2
              root3 = half * SQRT( three ) * ABS( x )
              root1 = - c0 / ( root2 * root2 + root3 * root3 )
              GO TO 900
            END IF
          END IF
        END IF

!  deflate cubic

        b0 = - c0 / root1
        IF ( ABS( root1 ** 3 ) <= ABS( c0 ) ) THEN
          b1 = root1 + three * c2
        ELSE
          b1 = ( b0 - three * c1 ) / root1
        END IF
        CALL ROOTS_quadratic( b0, b1, one, epsmch, nroots_q,                   &
                              root2, root3, debug )
        nroots = nroots + nroots_q


!  3. Use Viete's method

      ELSE IF ( method == 3 ) THEN
        w = a2 / ( three * a3 )
        p = ( a1 / ( three * a3 ) - w ** 2 ) ** 3
        q = - half * ( two * w ** 3 - ( a1 * w - a0 ) / a3 )
        d = p + q ** 2

!  three real roots

        IF ( d < zero ) THEN
          s = ACOS( MIN( one, MAX( - one, q / SQRT( - p ) ) ) )
          p = two * ( - p ) ** onesixth
          nroots = 3
          root1 = p * COS( onethird * ( s + two * pi ) ) - w
          root2 = p * COS( onethird * ( s + four * pi ) ) - w
          root3 = p * COS( onethird * ( s + six * pi ) ) - w

!  one real root

        ELSE
          d = SQRT( d ) ; u1 = q + d ; u2 = q - d
          nroots = 1
          root1 = SIGN( ABS( u1 ) ** onethird, u1 ) +                          &
                  SIGN( ABS( u2 ) ** onethird, u2 ) - w
        END IF

!  4. Compute the roots as the eigenvalues of the relevant compainion matrix

      ELSE
        H( 1, 1 ) = zero ; H( 2, 1 ) = one ; H( 3, 1 ) = zero
        H( 1, 2 ) = zero ; H( 2, 2 ) = zero ; H( 3, 2 ) = one
        H( 1, 3 ) = - a0 / a3 ; H( 2, 3 ) = - a1 / a3 ; H( 3, 3 ) = - a2 / a3
        CALL HSEQR( 'E', 'N', 3, 1, 3, H, 3, ER, EI, ZZ, 1, WORK, 33, info )
        IF ( info /= 0 ) THEN
          IF ( debug ) WRITE( out,                                             &
         &   "( ' ** error return ', I0, ' from HSEQR in ROOTS_cubic' )" ) info
          nroots = 0
          RETURN
        END IF

!  count and record the roots

        nroots = COUNT( ABS( EI ) <= epsmch )
        IF ( nroots == 1 ) THEN
          IF (  ABS( EI( 1 ) ) <= epsmch ) THEN
            root1 = ER( 1 )
          ELSE IF (  ABS( EI( 2 ) ) <= epsmch ) THEN
            root1 = ER( 2 )
          ELSE
            root1 = ER( 3 )
          END IF
        ELSE
          root1 = ER( 1 ) ;  root2 = ER( 2 ) ;  root3 = ER( 3 )
        END IF
      END IF

!  reorder the roots

  900 CONTINUE
      IF ( nroots == 3 ) THEN
        IF ( root1 > root2 ) THEN
          a = root2 ; root2 = root1 ; root1 = a
        END IF
        IF ( root2 > root3 ) THEN
          a = root3
          IF ( root1 > root3 ) THEN
            a = root1 ; root1 = root3
          END IF
          root3 = root2 ; root2 = a
        END IF
        IF ( debug ) WRITE( out, "( ' 3 real roots ' )" )
      ELSE IF ( nroots == 2 ) THEN
        IF ( debug ) WRITE( out, "( ' 2 real roots ' )" )
      ELSE
        IF ( debug ) WRITE( out, "( ' 1 real root ' )" )
      END IF

!  perfom a Newton iteration to ensure that the roots are accurate

      p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      pprime = ( three * a3 * root1 + two * a2 ) * root1 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime
        root1 = root1 - p / pprime
        p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 1, root1, p

      IF ( nroots == 3 ) THEN
        p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        pprime = ( three * a3 * root2 + two * a2 ) * root2 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime
          root2 = root2 - p / pprime
          p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 2, root2, p

        p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        pprime = ( three * a3 * root3 + two * a2 ) * root3 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 3, root3, p, - p / pprime
          root3 = root3 - p / pprime
          p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 3, root3, p
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4,         &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4 )


!  End of subroutine ROOTS_cubic

      END SUBROUTINE ROOTS_cubic

!-*-*-*-*-  G A L A H A D   T W O  _ N O R M   F U N C T I O N   -*-*-*-*-

       FUNCTION TWO_NORM( X )

!  Compute the l_2 norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: TWO_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         TWO_NORM = NRM2( n, X, 1 )
       ELSE
         TWO_NORM = zero
       END IF
       RETURN

!  End of function TWO_NORM

       END FUNCTION TWO_NORM

!-*-*-*-*-*-  End of R A L _ N L L S _ D T R S  double  M O D U L E  *-*-*-*-*-

   END MODULE RAL_NLLS_DTRS_double
