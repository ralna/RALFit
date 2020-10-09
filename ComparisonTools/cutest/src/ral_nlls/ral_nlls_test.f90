!  Dummy RAL_NLLS for testing ral_nlls_main interface to CUTEst
!  Nick Gould, 6th October 2015

module ral_nlls_workspaces
    
  implicit none 
  
  private
  
  ! define derived types and subroutines for workspace arrays.

  Integer, Parameter, Public          :: wp = kind(1.0d0)
  Real (Kind = wp), Parameter, Public :: epsmch = epsilon(1.0_wp)
  Real (Kind = wp), Parameter, Public :: toltm4= epsmch**(0.25_wp)
  Real (Kind = wp), Parameter, Public :: toltm8 = sqrt(epsmch)
  Real (Kind = wp), Parameter, Public :: toltm10 = sqrt(epsmch)/100.0_wp
  Real (Kind = wp), Parameter, Public :: toltm12= epsmch**(0.75_wp)
  Real (Kind = wp), Parameter, Public :: cmptol = toltm4/10.0_wp
  Real (Kind = wp), Parameter, Public :: infinity = HUGE(1.0_wp)
  INTEGER, PARAMETER, public          :: history_max = 100

  ! Error constants
  Integer, Parameter, Public :: NLLS_ERROR_MAXITS                   =   -1
  Integer, Parameter, Public :: NLLS_ERROR_EVALUATION               =   -2
  Integer, Parameter, Public :: NLLS_ERROR_UNSUPPORTED_MODEL        =   -3
  Integer, Parameter, Public :: NLLS_ERROR_FROM_EXTERNAL            =   -4
  Integer, Parameter, Public :: NLLS_ERROR_UNSUPPORTED_METHOD       =   -5
  Integer, Parameter, Public :: NLLS_ERROR_ALLOCATION               =   -6
  Integer, Parameter, Public :: NLLS_ERROR_MAX_TR_REDUCTIONS        =   -7
  Integer, Parameter, Public :: NLLS_ERROR_X_NO_PROGRESS            =   -8
  Integer, Parameter, Public :: NLLS_ERROR_BAD_TR_STRATEGY          =  -10
  Integer, Parameter, Public :: NLLS_ERROR_FIND_BETA                =  -11
  Integer, Parameter, Public :: NLLS_ERROR_BAD_SCALING              =  -12
  Integer, Parameter, Public :: NLLS_ERROR_WORKSPACE_ERROR          =  -13
  Integer, Parameter, Public :: NLLS_ERROR_UNSUPPORTED_TYPE_METHOD  =  -14
  Integer, Parameter, Public :: NLLS_ERROR_WRONG_INNER_METHOD       =  -15
  Integer, Parameter, Public :: NLLS_ERROR_INITIAL_GUESS            =  -16
  Integer, Parameter, Public :: NLLS_ERROR_UNSUPPORTED_LINESEARCH   =  -17
  Integer, Parameter, Public :: NLLS_ERROR_BAD_BOX_BOUNDS           =  -18


  ! dogleg errors
  Integer, Parameter, Public :: NLLS_ERROR_DOGLEG_MODEL             = -101
  ! AINT errors
  Integer, Parameter, Public :: NLLS_ERROR_AINT_EIG_IMAG            = -201
  Integer, Parameter, Public :: NLLS_ERROR_AINT_EIG_ODD             = -202
  ! More-Sorensen errors
  Integer, Parameter, Public :: NLLS_ERROR_MS_MAXITS                = -301
  Integer, Parameter, Public :: NLLS_ERROR_MS_TOO_MANY_SHIFTS       = -302
  Integer, Parameter, Public :: NLLS_ERROR_MS_NO_PROGRESS           = -303
  ! DTRS errors
  ! Tensor model errors
  Integer, Parameter, Public :: NLLS_ERROR_NO_SECOND_DERIVATIVES    = -401

  ! Linesearch errors
  Integer, Parameter, Public :: NLLS_ERROR_PG_STEP                  = -501

  ! Misc errors
  Integer, Parameter, Public :: NLLS_ERROR_PRINT_LEVEL              = -900
  Integer, Parameter, Public :: NLLS_ERROR_UNEXPECTED               = -999

    TYPE, PUBLIC :: NLLS_options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! P R I N T I N G   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   general output occurs on stream out

     INTEGER :: out = 6

!    The level of output required.
!    print_level = 0 gives no output (DEFAULT),
!    print_level = 1 gives a summary at the end of the solve,
!    print_level = 2 prints one-line summary for every iteration,
!    print_level = 3 same as 2 but add more details,
!    print_level = 4 same as 3 & also prints one-line inner iteration info,
!    print_level = 5 same as 4 with very verbose (debugging) output.

     INTEGER :: print_level = 0

!    Print all the options and their values
     Logical :: print_options = .False.

!    Print by default banner header every 30 iterations
     Integer :: print_header = 30

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

!   specify the model used. Possible values are
!
!      0  dynamic (*not yet implemented*)
!      1  Gauss-Newton (no 2nd derivatives)
!      2  second-order (exact Hessian)
!      3  hybrid (using Madsen, Nielsen and Tingleff's method)
!      4  tensor model (model must = 2)

     INTEGER :: model = 3

!   trust region or regularization?
!
!      1  trust region method
!      2  regularization

     INTEGER :: type_of_method = 1

!   specify the method used to solve the trust-region sub problem
!      1 Powell's dogleg
!      2 AINT method (of Yuji Nat.)
!      3 More-Sorensen
!      4 Galahad's DTRS

     INTEGER :: nlls_method = 4

!   allow the algorithm to use a different subproblem solver if one fails

     LOGICAL :: allow_fallback_method = .true.

     
!  which linear least squares solver should we use?

     INTEGER :: lls_solver = 1

!   overall convergence tolerances. The iteration will terminate when the
!   norm of the gradient of the objective function is smaller than
!      MAX( %stop_g_absolute, %stop_g_relative * norm of the initial gradient)
!   or if the norm of objective function is smaller than
!      MAX( %stop%f_absolute, %stop_f_relative * initial norm of the function)
!     or if the step is less than %stop_s

     REAL ( KIND = wp ) :: stop_g_absolute = 1.0e-5_wp
     REAL ( KIND = wp ) :: stop_g_relative = 1.0e-8_wp
     REAL ( KIND = wp ) :: stop_f_absolute = 1.0e-8_wp
     REAL ( KIND = wp ) :: stop_f_relative = 1.0e-8_wp
     REAL ( KIND = wp ) :: stop_s = epsmch


!   should we scale the initial trust region radius?

     integer :: relative_tr_radius = 0

!   if relative_tr_radius == 1, then pick a scaling parameter
!   Madsen, Nielsen and Tingleff say pick this to be 1e-6, say, if x_0 is good,
!   otherwise 1e-3 or even 1 would be good starts...

     real (kind = wp) :: initial_radius_scale = 1.0_wp

!   if relative_tr_radius /= 1, then set the
!   initial value for the trust-region radius (-ve => ||g_0||)

     REAL ( KIND = wp ) :: initial_radius = 100.0_wp

!   for the newton tensor model, allow a base tr raidius to allow an inherent
!   regularization in the problem that can't be changed
!   ( so we minimize 0.5 * (\sum (f_i)^2 + sigma_k ||s||^2) ), using another reg parameter
!   on top of this
!   (undocumented control variable)

     REAL ( KIND = wp ) :: base_regularization = 0.0_wp

     ! allow inherently the solution of a problem of the form
     !  min_x 1/2 ||r(x)||^2 + regularization_term * 1/ regularization_power * ||x||^regularization_term
     !
     ! this is done in two ways:
     !
     ! ** p = 2 **
     ! in this case, we solve a problem of the form
     !          min 0.5 || f(x) ||**2, where
     ! f:R^(n) -> R^(n+m)
     ! f_i(x) = r_i(x), i = 1,m
     ! f_i(x) = sqrt( regularization_term ) x_j, i = m + j, where j = 1,n
     ! This is implemented implicitly by updating
     !  ||f||**2 = ||r||**2 + regularization_term * ||x||**2
     !  J_f^Tf = J^Tr + regularization_term * x
     !  J_f^T J_f = J^T J + regularization_term * I
     !  md_f = md + 0.5 * regularization_term * ||x + d||**2
     !
     ! ** p .ne. 2 **
     ! here we solve a problem of the form
     !         min 0.5 || g(x) ||**2, where
     ! g:R^n -> R^(n+1)
     ! g_i(x) = r_i(x), i = 1,m
     ! g_i(x) = [(2*weight/power)**0.5 ] * ||x||**(power/2), i = m+1
     ! This is implemented implicitly by updating
     !  ||g||**2 = ||r||**2 + (2*weight/power) * ||x||**power
     !  J_g^T g = J^T r + weight ||x||**(power-2) x
     !  J_g^T J_g = J^T J + (weight * power / 2) * ||x||**(power-4) x x^T
     !  md_g = md + 0.5 * ||x||**(power-4) * weight *
     !             ( (2/power)**0.5 x^Tx + (power/2)**0.5 x^Td )**2
     !  and, if the full hessian was used
     !  md_g = md_g + weight * ||x||**(power-4)( x^Tx d^td + (d^tx)**2)
     integer :: regularization = 0
     REAL ( KIND = wp ) :: regularization_term = 1.0e-2_wp
     REAL ( KIND = wp ) :: regularization_power = 2.0_wp


!   maximum permitted trust-region radius

     REAL ( KIND = wp ) :: maximum_radius = 10.0_wp ** 8

!   a potential iterate will only be accepted if the actual decrease
!    f - f(x_new) is larger than %eta_successful times that predited
!    by a quadratic model of the decrease. The trust-region radius will be
!    increased if this relative decrease is greater than %eta_very_successful
!    but smaller than %eta_too_successful

     REAL ( KIND = wp ) :: eta_successful = 10.0_wp ** ( - 8 )! 10.0_wp ** ( - 8 )
     REAL ( KIND = wp ) :: eta_success_but_reduce = 10.0_wp ** ( - 8 ) !0.25_wp
     REAL ( KIND = wp ) :: eta_very_successful = 0.9_wp!0.75_wp!point9
     REAL ( KIND = wp ) :: eta_too_successful = 2.0_wp

!   on very successful iterations, the trust-region radius will be increased by
!    the factor %radius_increase, while if the iteration is unsuccessful, the
!    radius will be decreased by a factor %radius_reduce but no more than
!    %radius_reduce_max

     REAL ( KIND = wp ) :: radius_increase = 2.0_wp
     REAL ( KIND = wp ) :: radius_reduce = 0.5_wp
     REAL ( KIND = wp ) :: radius_reduce_max = 0.0625_wp

! Trust region update strategy
!    1 - usual step function
!    2 - continuous method of Hans Bruun Nielsen (IMM-REP-1999-05)
     integer :: tr_update_strategy = 1

!   if model=7, then the value with which we switch on second derivatives

     real ( kind = wp ) :: hybrid_switch = 0.1_wp

!   shall we use explicit second derivatives, or approximate using a secant
!   method

     LOGICAL :: exact_second_derivatives = .false.

!   use a factorization (dsyev) to find the smallest eigenvalue for the subproblem
!    solve? (alternative is an iterative method (dsyevx)
     LOGICAL :: subproblem_eig_fact = .FALSE. ! undocumented....

     ! use eigendecomposition in subproblem solve?
     LOGICAL :: use_ews_subproblem = .TRUE.

     ! This forces to call min_eig_symm without previously calling solve_spd_nocopy
     ! This option is used for code coverage and can be hidden from user.
     Logical :: force_min_eig_symm = .FALSE.

!   scale the variables?
!   0 - no scaling
!   1 - use the scaling in GSL (W s.t. W_ii = ||J(i,:)||_2^2)
!       tiny values get set to 1.0_wp
!   2 - scale using the approx to the Hessian (W s.t. W = ||H(i,:)||_2^2
     INTEGER :: scale = 1
     REAL(wp) :: scale_max = 1.0e11_wp
     REAL(wp) :: scale_min = 1.0e-11_wp
     LOGICAL :: scale_trim_min = .TRUE.
     LOGICAL :: scale_trim_max = .TRUE.
     LOGICAL :: scale_require_increase = .FALSE.

     logical :: setup_workspaces = .true.
     logical :: remove_workspaces = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! M O R E - S O R E N S E N   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer  :: more_sorensen_maxits = 500
    real(wp) :: more_sorensen_shift = 1.0e-13_wp
    real(wp) :: more_sorensen_tiny = 10.0_wp * epsmch
    real(wp) :: more_sorensen_tol = 1.0e-3_wp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! H Y B R I D   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what's the tolerance such that ||J^T f || < tol * 0.5 ||f||^2 triggers a switch
    real(wp) :: hybrid_tol = 2.0_wp!0.02

! how many successive iterations does the above condition need to hold before we switch?
    integer  :: hybrid_switch_its = 1!3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! T E N S O R   M O D E L   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what regularization should we use?
    real(wp) :: reg_order = -1.0_wp

! which method shall we use to solve the inner problem?
! 1 - add in a base regularization parameter
! 2 - minimize a modified cost functional
    integer :: inner_method = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! O U T P U T   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Shall we output progess vectors at termination of the routine?
     logical :: output_progress_vectors = .false.

     logical :: update_lower_order = .true.

     ! C or Fortran storage of Jacobian?
     logical :: Fortran_Jacobian = .true.
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! B O X   B O U N D   O P T I O N S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Memory size for the non-monotone linesearch
     Integer :: box_nFref_max = 1
!    Kanzow sufficient decrease ratio (eq 25) Kanzow 2004
     Real(Kind=wp) :: box_gamma = 0.99995_wp
     Real(Kind=wp) :: box_decmin = epsmch
!    Magic number to consider box bound (+/-) infinity
     Real(Kind=wp) :: box_bigbnd = 1.0e20_wp
!    Wolfe descent condition (0<\sigma1<1/2), curvature condition (0<\sigma2)
     Real(Kind=wp) :: box_wolfe_descent = 1.0e-4_wp
     Real(Kind=wp) :: box_wolfe_curvature = 0.9_wp
!    Tolerance to consider projected dir a descent direction
!    See LS STEP Section 4 p392 Kanzow 2014
     Real(Kind=wp) :: box_kanzow_power = 2.1_wp
!    sqrt(mcheps)
     Real(Kind=wp) :: box_kanzow_descent = toltm8 
!    sqrt(mcheps)
     Real(Kind=wp) :: box_quad_model_descent = toltm8
!    Take projected TR step when TR test is Ok?
!    True  => take step
!    False => force a LS or PG step
     Logical       :: box_tr_test_step = .True.
!    Take projected  TR step when Wolfe test is Ok?
!    True  => take step
!    False => force a LS or PG step
     Logical       :: box_wolfe_test_step = .True.
!    Threshold to determine if the projection of TR direction
!    is too severe 0<tau_min<1
     Real(Kind=wp) :: box_tau_min = 0.1_wp
!    tau >= tau_descent in order to test for descent
     Real(Kind=wp) :: box_tau_descent = 1.0e-5_wp
!    Max times TR iterations can fail without passing the various
!    descent tests: 2? 3? 5? Ignored when proj(x)==x
     Integer       :: box_max_ntrfail = 2
!    Number of consecutive times quadratic model matches f(x_k+1)
!    required before setting initial alpha step for PG step equal
!    to scale_alpha*alpha_k-1
     Integer       :: box_quad_match = 2
!    Initial step scale (if quad_i >= box_quad_i)
     Real(Kind=wp) :: box_alpha_scale = 1.0_wp
!    Scaling factor to use when updating Delta from LS/PG step
     Real(Kind=wp) :: box_Delta_scale = 2.0_wp
     Real(Kind=wp) :: box_tau_wolfe = 0.3_wp
     Real(Kind=wp) :: box_tau_tr_step = 0.3_wp
     Integer       :: box_ls_step_maxit = 20
!    LS type: 1 => Dennis-Schnabel; 2 => Hager-Zhang
     Integer       :: box_linesearch_type = 1
  END TYPE nlls_options

!  - - - - - - - - - - - - - - - - - - - - - - -
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

  TYPE, PUBLIC :: nlls_inform

!  return status
!  (see ERROR type for descriptions)
     INTEGER :: status = 0

! error message

     CHARACTER ( LEN = 80 ) :: error_message = REPEAT( ' ', 80 )


!  the status of the last attempted allocation/deallocation

     INTEGER :: alloc_status = 0

!  the name of the array for which an allocation/deallocation error ocurred

     CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  the total number of iterations performed

     INTEGER :: iter = 0

!  the number of inner iterations performed

     INTEGER :: inner_iter = 0

!  exit status of last inner iteration process

     LOGICAL :: inner_iter_success = .False.

!  the total number of CG iterations performed

!$$     INTEGER :: cg_iter = 0

!  the total number of evaluations of the objective function

     INTEGER :: f_eval = 0

!  the total number of evaluations of the gradient of the objective function

     INTEGER :: g_eval = 0

!  the total number of evaluations of the Hessian of the objective function
!  using eval_hf
     
     INTEGER :: h_eval = 0

!  the total number of evaluations of the Hessian of the objective function
!  using eval_hp

     INTEGER :: hp_eval = 0

!  test on the size of f satisfied?

     integer :: convergence_normf = 0

!  test on the size of the gradient satisfied?

     integer :: convergence_normg = 0

!  test on the size of the step satisfied?

     integer :: convergence_norms = 0

!  vector of residuals

     real(wp), allocatable :: resvec(:)

!  vector of gradients

     real(wp), allocatable :: gradvec(:)

!  the value of the objective function at the best estimate of the solution
!   determined by NLLS_solve

     REAL ( KIND = wp ) :: obj = infinity

!  the norm of the gradient of the objective function at the best estimate
!   of the solution determined by NLLS_solve

     REAL ( KIND = wp ) :: norm_g = infinity

! the norm of the gradient, scaled by the norm of the residual

     REAL ( KIND = wp ) :: scaled_g = infinity

! error returns from external subroutines

     INTEGER :: external_return = 0

! name of external program that threw and error

     CHARACTER ( LEN = 80 ) :: external_name = REPEAT( ' ', 80 )

! Step size of last iteration (w%normd)
     Real(Kind=wp) :: step = infinity

! LS Step iterations
     Integer :: ls_step_iter = 0
     Integer :: f_eval_ls = 0
     Integer :: g_eval_ls = 0
! PG Step iterations
     Integer :: pg_step_iter = 0
     Integer :: f_eval_pg = 0
     Integer :: g_eval_pg = 0

     
  END TYPE nlls_inform

  type, public :: params_base_type
     ! deliberately empty
  end type params_base_type

  abstract interface
     subroutine eval_hf_type(status, n, m, x, f, h, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: f
       double precision, dimension(*), intent(out) :: h
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hf_type
  end interface

  abstract interface
     subroutine eval_hp_type(status, n, m, x, y, hp, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: y
       double precision, dimension(*), intent(out) :: hp
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hp_type
  end interface

  type, extends( params_base_type ), public :: tensor_params_type
     real(wp), dimension(:), allocatable :: f
     real(wp), dimension(:), allocatable :: x
     real(wp), dimension(:), allocatable :: J
     real(wp), dimension(:,:,:), allocatable :: Hi
     real(wp) :: Delta
     real(wp) :: p
     integer :: m
     integer :: m1 = 0
     integer :: extra = 0
     procedure( eval_hf_type ), pointer, nopass :: eval_HF
     procedure( eval_hp_type ), pointer, nopass :: eval_HP
     logical :: eval_hp_provided = .false.
     class( params_base_type ), pointer :: parent_params
     type( tenJ_type ), pointer :: tenJ
  end type tensor_params_type

    type, public :: box_type
     ! Does the problem have a box?
     logical :: has_box = .false.
     Real(Kind=wp), Allocatable :: blx(:), bux(:), pdir(:), normFref(:), sk(:), g(:)
     ! projection changed the direction? d /= P(d)?
     Logical :: prjchd = .false.
     ! Convergence metrics
     Real(Kind=wp) :: normPD, gtd
     ! Memory for HZLS (LS)
     Real(Kind=wp) :: sksk, skyk, quad_c, quad_q, normFold
     ! Consecutive times quadratic model is accurate
     Integer       :: quad_i = 0
     ! Memory for nonmonotone LS
     Integer       :: nFref = 0
  end type box_type
  
    type, public :: max_eig_work ! workspace for subroutine max_eig
     logical :: allocated = .false.
     real(wp), allocatable :: alphaR(:), alphaI(:), beta(:), vr(:,:)
     real(wp), allocatable :: work(:), ew_array(:)
     integer, allocatable :: nullindex(:)
     logical, allocatable :: vecisreal(:)
     integer :: nullevs_cols
     real(wp), allocatable :: nullevs(:,:)
  end type max_eig_work

  type, public :: minus_solve_general_work ! workspace for subroutine solve_general
     logical :: allocated = .false.
     real(wp), allocatable :: A(:,:)
     integer, allocatable :: ipiv(:)
  end type minus_solve_general_work

  type, public :: evaluate_model_work ! workspace for subroutine evaluate_model
     logical :: allocated = .false.
     real(wp), allocatable :: Jd(:), dH(:), Hd(:), dHd(:)
  end type evaluate_model_work

  type, public :: solve_LLS_work ! workspace for subroutine solve_LLS
     logical :: allocated = .false.
     real(wp), allocatable :: temp(:), work(:), Jlls(:)
  end type solve_LLS_work

  type, public :: min_eig_symm_work ! workspace for subroutine min_eig_work
     logical :: allocated = .false.
     real(wp), allocatable :: A(:,:), work(:), ew(:)
     integer, allocatable :: iwork(:), ifail(:)
  end type min_eig_symm_work

  type, public :: all_eig_symm_work ! workspace for subroutine all_eig_symm
     logical :: allocated = .false.
     real(wp), allocatable :: work(:)
     ! This will only be allocated with ew(n) to comply with LAPACK interface
     real(wp), allocatable :: ew(:)
  end type all_eig_symm_work

  type, public :: generate_scaling_work ! workspace for subrouine generate_scaling
     logical :: allocated = .false.
     real(wp), allocatable :: diag(:)
     real(wp), allocatable :: ev(:,:)
     real(wp), allocatable :: tempvec(:)
     type( all_eig_symm_work ) :: all_eig_symm_ws
  end type generate_scaling_work

  type, public :: solve_newton_tensor_work ! a workspace for solve_newton_tensor
     logical :: allocated = .false.
     real(wp), allocatable :: model_tensor(:)
!     type( generate_scaling_ws ) ::
     type( tensor_params_type ) :: tparams
     type( nlls_options ) :: tensor_options
     integer :: m_in
  end type solve_newton_tensor_work

  type, public :: solve_galahad_work ! workspace for subroutine dtrs_work
     logical :: allocated = .false.
     real(wp), allocatable ::ev(:,:), ew(:), v_trans(:), d_trans(:), scale_c(:),&
       scale_h(:)
     type( all_eig_symm_work ) :: all_eig_symm_ws
  end type solve_galahad_work

  type, public :: regularization_solver_work ! workspace for subroutine regularization_solver
     logical :: allocated = .false.
     real(wp), allocatable :: AplusSigma(:,:),LtL(:,:)
  end type regularization_solver_work

  type, public :: more_sorensen_work ! workspace for subroutine more_sorensen
     logical :: allocated = .false.
     real(wp), allocatable :: LtL(:,:), AplusSigma(:,:)
     real(wp), allocatable :: q(:), y1(:)
     type( min_eig_symm_work ) :: min_eig_symm_ws
     real(wp), allocatable :: norm_work(:)
  end type more_sorensen_work

  type, public :: AINT_tr_work ! workspace for subroutine AINT_tr
     logical :: allocated = .false.
     type( max_eig_work ) :: max_eig_ws
     type( evaluate_model_work ) :: evaluate_model_ws
     type( minus_solve_general_work ) :: minus_solve_general_ws
     !       type( solve_spd_work ) :: minus_solve_spd_ws
     REAL(wp), allocatable :: LtL(:,:), B(:,:), p0(:), p1(:)
     REAL(wp), allocatable :: M0(:,:), M1(:,:), y(:), gtg(:,:), q(:)
     REAL(wp), allocatable :: M0_small(:,:), M1_small(:,:)
     REAL(wp), allocatable :: By_hardcase(:,:)
  end type AINT_tr_work

  type, public :: dogleg_work ! workspace for subroutine dogleg
     logical :: allocated = .false.
     type( solve_LLS_work ) :: solve_LLS_ws
     type( evaluate_model_work ) :: evaluate_model_ws
     real(wp), allocatable :: d_sd(:), d_gn(:), ghat(:), Jg(:)
  end type dogleg_work

  type, public :: calculate_step_work ! workspace for subroutine calculate_step
     logical :: allocated = .false.
     real(wp), allocatable :: A(:,:), xxt(:,:)
     real(wp), allocatable :: v(:), scale(:), extra_scale(:)
     real(wp) :: reg_order = 2.0_wp ! reg. by + 1/p || \sigma || ** p
     type( AINT_tr_work ) :: AINT_tr_ws
     type( dogleg_work ) :: dogleg_ws
     type( solve_newton_tensor_work ) :: solve_newton_tensor_ws
     type( more_sorensen_work ) :: more_sorensen_ws
     type( solve_galahad_work ) :: solve_galahad_ws
     type( regularization_solver_work ) :: regularization_solver_ws
     type( evaluate_model_work) :: evaluate_model_ws
     type( generate_scaling_work ) :: generate_scaling_ws
  end type calculate_step_work

  type, public :: tenJ_type ! workspace for evaltensor_J
     logical :: allocated = .false.
     real(wp), allocatable :: Hs(:,:), Js(:) ! work arrays for evaltensor_f
     real(wp), allocatable :: H(:,:)       ! and another....
     real(wp), allocatable :: stHs(:)      ! yet another....
  end type tenJ_type

  type, public :: NLLS_workspace ! all workspaces called from the top level
     logical :: allocated = .false.
     integer :: first_call = 1
     integer :: iter = 0
     real(wp) :: normF0, normJF0, normF, normJF
     real(wp) :: normJFold, normJF_Newton
     real(wp) :: Delta
     real(wp) :: norm_2_d ! 2-norm of d
     real(wp) :: norm_S_d ! S-norm of d, where S is the scaling
     logical :: use_second_derivatives = .false.
     integer :: hybrid_count = 0
     real(wp) :: hybrid_tol = 1.0_wp
     real(wp), allocatable :: fNewton(:), JNewton(:), XNewton(:)
     real(wp), allocatable :: J(:)
     real(wp), allocatable :: f(:), fnew(:)
     real(wp), allocatable :: hf(:), hf_temp(:)
     real(wp), allocatable :: d(:), g(:), Xnew(:)
     real(wp), allocatable :: y(:), y_sharp(:), g_old(:), g_mixed(:)
     real(wp), allocatable :: ysharpSks(:), Sks(:)
     real(wp), allocatable :: resvec(:), gradvec(:)
     real(wp), allocatable :: Wf(:)
     type ( calculate_step_work ) :: calculate_step_ws
     type ( box_type ) :: box_ws
     real(wp) :: tr_nu = 2.0_wp
     integer :: tr_p = 3
     type (tenJ_type ) :: tenJ
     type (NLLS_workspace), Pointer :: iw_ptr => NULL()
  end type NLLS_workspace

  
end module ral_nlls_workspaces

module ral_nlls_internal 
     
  use iso_c_binding
  use ral_nlls_workspaces

  implicit none

  private
  
  INTEGER, PARAMETER :: error_dimensions = - 1
  INTEGER, PARAMETER :: error_workspace = - 2
  INTEGER, PARAMETER :: error_eval_F = - 3
  INTEGER, PARAMETER :: error_eval_J = - 4
  INTEGER, PARAMETER :: error_eval_HF = - 5
  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: f
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hf_type
  end interface

  abstract interface
     subroutine eval_hp_type(status, n, m, x, y, hp, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: y
       double precision, dimension(*), intent(out) :: hp
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hp_type     
  end interface

  public :: nlls_solve
  
contains

          SUBROUTINE NLLS_SOLVE( n, m, X, eval_F, eval_J, eval_HF,                    &
               params, options, inform, weights, eval_HP, &
               lower_bounds, upper_bounds)
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

     INTEGER( c_int ), INTENT( IN ) :: n, m
     REAL( c_double ), DIMENSION( n ), INTENT( INOUT ) :: X
     TYPE( Nlls_inform ), INTENT( OUT ) :: inform
     TYPE( Nlls_options ), INTENT( IN ) :: options
     procedure( eval_f_type ) :: eval_F
     procedure( eval_j_type ) :: eval_J
     procedure( eval_hf_type ) :: eval_HF
     class( params_base_type ), intent(inout) :: params
     real( wp ), dimension(m), intent(in), optional :: weights
     procedure( eval_hp_type ), optional :: eval_HP
     real( wp ), dimension( n ), intent(inout), optional :: lower_bounds
     real( wp ), dimension( n ), intent(inout), optional :: upper_bounds

     type ( NLLS_workspace ) :: w
     Type ( NLLS_workspace ), target :: inner_workspace
     
!  Interface blocks

!  Local variables

     INTEGER :: status, start_f, end_f, start_j, start_h, w_end
     INTEGER :: len_work_int, len_work_real
     INTEGER( c_int ), allocatable :: Work_int( : )
     REAL( c_double ), allocatable :: Work_real( : ) 
     
!  check input dimensions

     IF ( m <= 0 .OR. n <= 0 ) THEN
       status = error_dimensions
       GO TO 990
     END IF

!  partition the workspace
     allocate(Work_int(10))
     allocate(Work_real(m + n*(n + m)))
     
     start_f = 1
     start_j = start_f + m
     end_f = start_j - 1
     start_h = start_j + n * m
     w_end = start_h + n * n - 1

     IF ( w_end < len_work_real ) THEN
       status = error_workspace
       GO TO 990
     END IF

!  evaluate F

     CALL eval_F( status, n, m, X, WORK_real( start_f ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_F
       GO TO 990
     END IF
     inform%obj = 0.5_c_double * DOT_PRODUCT( WORK_real( start_f : end_f ),    &
                                              WORK_real( start_f : end_f ) )

!  evaluate J

     CALL eval_J( status, n, m, X, WORK_real( start_j ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_J
       GO TO 990
     END IF

!  evaluate HF

     CALL eval_HF( status, n, m, X, WORK_real( start_f ), WORK_real( start_h ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_HF
       GO TO 990
     END IF

 990 CONTINUE
     inform%status = status
     RETURN
   END SUBROUTINE NLLS_SOLVE
     
end module ral_nlls_internal


   MODULE RAL_NLLS_DOUBLE

     USE ISO_C_BINDING
     use ral_nlls_internal
     use ral_nlls_workspaces, only : params_base_type, nlls_options, &  
                                     nlls_inform, nlls_workspace
     
     IMPLICIT none

     private

!!$     INTEGER, PARAMETER :: wp = KIND( 1.0d0 )
!!$     INTEGER, PARAMETER :: long = SELECTED_INT_KIND( 8 )
!!$

!!$     
!!$     real (kind = wp), parameter :: tenm3 = 1.0e-3
!!$     real (kind = wp), parameter :: tenm5 = 1.0e-5
!!$     real (kind = wp), parameter :: tenm8 = 1.0e-8
!!$     real (kind = wp), parameter :: epsmch = epsilon(1.0_wp)
!!$     real (kind = wp), parameter :: hundred = 100.0
!!$     real (kind = wp), parameter :: ten = 10.0
!!$     real (kind = wp), parameter :: point9 = 0.9
!!$     real (kind = wp), parameter :: zero = 0.0
!!$     real (kind = wp), parameter :: one = 1.0
!!$     real (kind = wp), parameter :: two = 2.0
!!$     real (kind = wp), parameter :: half = 0.5
!!$     real (kind = wp), parameter :: sixteenth = 0.0625
  
     
  ABSTRACT INTERFACE
     SUBROUTINE eval_F_type( status, n, m, X, F , params )
       import :: params_base_type
       implicit none
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n, m
       double precision, DIMENSION( * ), INTENT( IN ) :: X
       double precision, DIMENSION( * ), INTENT( OUT ) :: F
       class( params_base_type ), intent( inout ) :: params
     END SUBROUTINE eval_F_type
  END INTERFACE

  ABSTRACT INTERFACE
     SUBROUTINE eval_j_type( status, n, m, X, J, params )
       USE ISO_C_BINDING
       import :: params_base_type
       INTEGER ( c_int ), INTENT( OUT ) :: status
       INTEGER ( c_int ), INTENT( IN ) :: n, m
       REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
       REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: J
       class( params_base_type ), intent( inout ) :: params
     END SUBROUTINE eval_j_type
  END INTERFACE

  ABSTRACT INTERFACE
     SUBROUTINE eval_HF_type( status, n, m, X, F, H, params ) 
       USE ISO_C_BINDING
       import :: params_base_type
       INTEGER ( c_int ), INTENT( OUT ) :: status
       INTEGER ( c_int ), INTENT( IN ) :: n, m
       REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
       REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: F
       REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: H
       class( params_base_type ), intent( inout ) :: params
     END SUBROUTINE eval_HF_type
  END INTERFACE

  public :: nlls_solve
  public :: nlls_options, nlls_inform, nlls_workspace
  public :: params_base_type

CONTAINS

  ! empty

   END MODULE RAL_NLLS_DOUBLE


