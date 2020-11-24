! ral_nlls_workspaces :: module to keep all the workspaces

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

     ! This forces to call min_eig_symm without previously calling minus_solve_spd_nocopy
     ! This option is used for code coverage and can be hidden from user.
     Logical :: force_min_eig_symm = .FALSE.

!   scale the variables?
!   0 - no scaling
!   1 - use the scaling in GSL (W s.t. W_ii = ||J(i,:)||_2^2)
!       tiny values get set to one
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
!       Save covariance matrix type
!         0: None
!         1: C = sigma^2 * inv(J^T J)
!         2: only diagonal of C
!         3: only J^T J
     Integer       :: save_covm = 0

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

!  the norm of the gradient or projected gradient of the objective function 
!  at the current best estimate of the solution determined by NLLS_solve

     REAL ( KIND = wp ) :: norm_g = infinity

! the norm of the gradient (or projected gradient), scaled by the norm of 
! the residual

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

! COVARIANCE MATRIX (VARIANCE vector or J^T*J)
     real(wp), allocatable :: cov(:,:)
     real(wp), allocatable :: var(:)

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
     logical :: Fortran_Jacobian
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

  type, public :: minus_solve_general_work ! workspace for subroutine minus_solve_general
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
     !       type( minus_solve_spd_work ) :: minus_solve_spd_ws
     REAL(wp), allocatable :: LtL(:,:), B(:,:), p0(:), p1(:)
     REAL(wp), allocatable :: M0(:,:), M1(:,:), y(:), gtg(:,:), q(:)
     REAL(wp), allocatable :: M0_small(:,:), M1_small(:,:)
     REAL(wp), allocatable :: y_hardcase(:,:)
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

  public :: setup_workspaces, remove_workspaces
  public :: setup_workspace_dogleg, setup_workspace_AINT_tr
  public :: setup_workspace_more_sorensen, setup_workspace_solve_galahad
  public :: setup_workspace_regularization_solver
  public :: setup_bounds_type, remove_workspace_bounds

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                       !!
  !! W O R K S P A C E   S E T U P   S U B R O U T I N E S !!
  !!                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine setup_workspaces(workspace,n,m,options,inform)
  implicit none
    type( NLLS_workspace ), intent(inout) :: workspace
    type( nlls_options ), intent(in) :: options
    integer, intent(in) :: n,m
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    if (.not. allocated(workspace%y)) then
       allocate(workspace%y(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
       workspace%y = 0.0_wp
    end if

    if (.not. allocated(workspace%y_sharp)) then
       allocate(workspace%y_sharp(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
       workspace%y_sharp = 0.0_wp
    end if

    if (.not. options%exact_second_derivatives) then
       if (.not. allocated(workspace%g_old)) then
          allocate(workspace%g_old(n), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
       if (.not. allocated(workspace%g_mixed)) then
          allocate(workspace%g_mixed(n), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
       if (.not. allocated(workspace%Sks)) then
          allocate(workspace%Sks(n), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
       if (.not. allocated(workspace%ysharpSks)) then
          allocate(workspace%ysharpSks(n), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
    end if

    if( options%output_progress_vectors ) then
       if (.not. allocated(workspace%resvec)) then
          allocate(workspace%resvec(options%maxit+1), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
       if (.not. allocated(workspace%gradvec)) then
          allocate(workspace%gradvec(options%maxit+1), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
    end if

    if( .not. allocated(workspace%J)) then
       allocate(workspace%J(n*m), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%f)) then
       allocate(workspace%f(m), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%Wf)) then
       allocate(workspace%Wf(m), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%fnew)) then
       allocate(workspace%fnew(m), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%hf)) then
       allocate(workspace%hf(n*n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( options%model == 3 ) then
       if( .not. allocated(workspace%hf_temp)) then
          allocate(workspace%hf_temp(n*n), stat = inform%alloc_status)
          If (inform%alloc_status /= 0) Then
            inform%bad_alloc = 'setup_workspaces'
            inform%status = NLLS_ERROR_ALLOCATION
            goto 100
          End If
       end if
    end if

    if( .not. allocated(workspace%d)) then
       allocate(workspace%d(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%g)) then
       allocate(workspace%g(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if

    if( .not. allocated(workspace%Xnew)) then
       allocate(workspace%Xnew(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%bad_alloc = 'setup_workspaces'
         inform%status = NLLS_ERROR_ALLOCATION
         goto 100
       End If
    end if
    call setup_workspace_calculate_step(n,m,workspace%calculate_step_ws, &
         options, inform, workspace%tenJ, workspace%iw_ptr)
    if (inform%status/=0) goto 100

    workspace%allocated = .true.

100 continue

!   Clean-up if error encountered
    If (inform%alloc_status /= 0) Call remove_workspaces(workspace, options)

  end subroutine setup_workspaces

  recursive subroutine remove_workspaces(workspace,options)
    implicit none
    type( NLLS_workspace ), intent(inout) :: workspace
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated(workspace%y)) deallocate(workspace%y, stat=ierr_dummy)
    if(allocated(workspace%y_sharp)) deallocate(workspace%y_sharp, stat=ierr_dummy)
    if(allocated(workspace%g_old)) deallocate(workspace%g_old, stat=ierr_dummy)
    if(allocated(workspace%g_mixed)) deallocate(workspace%g_mixed, stat=ierr_dummy)
    if(allocated(workspace%Sks)) deallocate(workspace%Sks, stat=ierr_dummy)
    if(allocated(workspace%ysharpSks)) deallocate(workspace%ysharpSks, stat=ierr_dummy)

    if(allocated(workspace%resvec)) deallocate(workspace%resvec, stat=ierr_dummy)
    if(allocated(workspace%gradvec)) deallocate(workspace%gradvec, stat=ierr_dummy)

    if(allocated(workspace%fNewton)) deallocate(workspace%fNewton, stat=ierr_dummy )
    if(allocated(workspace%JNewton)) deallocate(workspace%JNewton, stat=ierr_dummy )
    if(allocated(workspace%XNewton)) deallocate(workspace%XNewton, stat=ierr_dummy )

    if(allocated(workspace%J)) deallocate(workspace%J, stat=ierr_dummy )
    if(allocated(workspace%f)) deallocate(workspace%f, stat=ierr_dummy )
    if(allocated(workspace%Wf)) deallocate(workspace%Wf, stat=ierr_dummy )
    if(allocated(workspace%fnew)) deallocate(workspace%fnew, stat=ierr_dummy )
    if(allocated(workspace%hf)) deallocate(workspace%hf, stat=ierr_dummy )
    if(allocated(workspace%hf_temp)) deallocate(workspace%hf_temp, stat=ierr_dummy)
    if(allocated(workspace%d)) deallocate(workspace%d, stat=ierr_dummy )
    if(allocated(workspace%g)) deallocate(workspace%g, stat=ierr_dummy )
    if(allocated(workspace%Xnew)) deallocate(workspace%Xnew, stat=ierr_dummy )

    call remove_workspace_calculate_step(workspace%calculate_step_ws,&
         options,workspace%tenJ, workspace%iw_ptr)

    call remove_workspace_bounds(workspace%box_ws)

    workspace%allocated = .false.

  end subroutine remove_workspaces


  subroutine setup_workspace_solve_newton_tensor(n,m,w,options,inform,tenJ,inner_workspace)
    implicit none
    integer, intent(in) :: n, m
    type( solve_newton_tensor_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform
    type( tenJ_type ), intent(InOut) :: tenJ
    type( NLLS_workspace ), intent(InOut) :: inner_workspace
    integer :: ierr_dummy

    inform%status = 0
    allocate(w%model_tensor(m),w%tparams%f(m),w%tparams%x(n),w%tparams%J(n*m), &
      w%tparams%Hi(n,n,m), tenJ%Hs(n,m), tenJ%Js(m),tenJ%H(n,n),tenJ%stHs(m),  &
      stat=inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%bad_alloc = 'setup_workspace_solve_newton_tensor'
      inform%status = NLLS_ERROR_ALLOCATION
      goto 100
    End If

    w%tparams%m = m

    ! copy options from those input
    w%tensor_options = options
    ! use a hybrid method for the inner loop
    w%tensor_options%model = 3
    w%tensor_options%maxit = 100
    w%tensor_options%reg_order = -1.0_wp
    w%tensor_options%output_progress_vectors = .false.
    select case (options%inner_method)
    case (1)
       w%tensor_options%type_of_method = 2
       w%tensor_options%nlls_method = 4
       !         w%tensor_options%radius_increase = 2.0_wp
       !         w%tensor_options%radius_reduce = 0.5_wp
       ! we seem to get better performance using
       ! straight Newton here (Why?)
       w%tensor_options%model = 2
       w%m_in = m
    case (2,5)
       w%tensor_options%model = 3
       w%tensor_options%type_of_method = 2 ! make this changable by the user
       w%tensor_options%nlls_method = 4
       !         w%tensor_options%radius_increase = 2.0_wp
       !         w%tensor_options%radius_reduce = 0.5_wp
       w%tensor_options%stop_g_absolute = 1.0e-10_wp
       w%tensor_options%stop_g_relative = 1.0e-10_wp
       w%tparams%m1 = m
       if (options%inner_method == 2) then
          w%m_in = m + 1
       else
          w%m_in = m + n
       end if
    case (3,4)
       w%tensor_options%model = 3
       w%tensor_options%type_of_method = 2
       w%tensor_options%nlls_method = 4
       !         w%tensor_options%stop_g_absolute = 1e-14
       !         w%tensor_options%stop_g_relative = 1e-14
       !          w%tensor_options%radius_increase = 2.0_wp
       !          w%tensor_options%radius_reduce = 0.5_wp
       w%tparams%m1 = m
       w%m_in = m
       if (options%inner_method == 3) then
          w%tensor_options%regularization = 1
          w%tensor_options%regularization_power = 2.0_wp
       else
          w%tensor_options%regularization = 2
          w%tensor_options%regularization_power = 3.0_wp
       end if
    case Default
       inform%status = NLLS_ERROR_WRONG_INNER_METHOD
       Go To 100
    end select

    ! setup/remove workspaces manually....
    w%tensor_options%remove_workspaces = .false.

    w%tensor_options%setup_workspaces = .false.
    call setup_workspaces(inner_workspace, n, w%m_in, w%tensor_options, inform)
    if (inform%status /= 0) GoTo 100

    w%allocated = .true.

100 continue
    if (inform%status /= 0) then
      If (Allocated(w%tparams%x)) Deallocate(w%tparams%x, stat=ierr_dummy)
      Call remove_workspace_solve_newton_tensor(w,options,tenJ,inner_workspace)
    end if
  end subroutine setup_workspace_solve_newton_tensor

  subroutine remove_workspace_solve_newton_tensor(w,options,tenJ,inner_workspace)
    implicit none
    type( solve_newton_tensor_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( tenJ_type ), Intent(InOut) :: tenJ
    type( NLLS_workspace ), Intent(InOut) :: inner_workspace
    Integer :: ierr_dummy

    if(allocated(w%model_tensor)) deallocate(w%model_tensor,stat=ierr_dummy)
    if(allocated(w%tparams%f)) deallocate(w%tparams%f,stat=ierr_dummy)
    if(allocated(w%tparams%J)) deallocate(w%tparams%J,stat=ierr_dummy)
    if(allocated(w%tparams%Hi)) deallocate(w%tparams%Hi,stat=ierr_dummy)
    if(allocated(tenJ%Hs)) deallocate(tenJ%Hs,stat=ierr_dummy)
    if(allocated(tenJ%Js)) deallocate(tenJ%Js,stat=ierr_dummy)
    if(allocated(tenJ%H)) deallocate(tenJ%H,stat=ierr_dummy)
    if(allocated(tenJ%stHs)) deallocate(tenJ%stHs,stat=ierr_dummy)
    call remove_workspaces(inner_workspace, w%tensor_options)

    w%allocated = .false.
  end subroutine remove_workspace_solve_newton_tensor

  recursive subroutine setup_workspace_calculate_step(n,m,w,options,inform,tenJ,inner_workspace)
    implicit none
    integer, intent(in) :: n, m
    type( calculate_step_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform
    type( tenJ_type ), intent(InOut) :: tenJ
    type( NLLS_workspace ), intent(InOut) :: inner_workspace

    inform%status = 0
    allocate(w%A(n,n), w%v(n),w%extra_scale(n),w%scale(n),stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%bad_alloc = 'setup_workspace_calculate_step'
      inform%status = NLLS_ERROR_ALLOCATION
      goto 100
    End If

    w%scale(:) = 0.0_wp

    call setup_workspace_evaluate_model(n,m,&
         w%evaluate_model_ws,options,inform)
    if (inform%status /= 0) goto 100
    if ( options%model == 4 ) then
       call setup_workspace_solve_newton_tensor(n,m,&
            w%solve_newton_tensor_ws,&
            options, inform, tenJ, inner_workspace)
       if (inform%status /= 0) goto 100
    else
       if ( options%type_of_method == 1) then
          select case (options%nlls_method)
          case (1) ! use the dogleg method
             call setup_workspace_dogleg(n,m,w%dogleg_ws, &
               options, inform)
             if (inform%status /= 0) goto 100
          case(2) ! use the AINT method
             call setup_workspace_AINT_tr(n,m,w%AINT_tr_ws, &
                  options, inform)
             if (inform%status /= 0) goto 100
          case(3) ! More-Sorensen
             call setup_workspace_more_sorensen(n,m,&
                  w%more_sorensen_ws,options,inform)
             if (inform%status /= 0) goto 100
          case (4) ! dtrs (Galahad)
             call setup_workspace_solve_galahad(n,m, &
                  w%solve_galahad_ws, options, inform)
             if (inform%status /= 0) goto 100
          Case Default
             inform%status = NLLS_ERROR_UNSUPPORTED_METHOD
             GoTo 100
          end select
       elseif (options%type_of_method == 2) then
          select case (options%nlls_method)
          case (3)
             call setup_workspace_regularization_solver(n,m, &
                  w%regularization_solver_ws, options, inform)
             if (inform%status /= 0) goto 100
          case (4)
             call setup_workspace_solve_galahad(n,m, &
                  w%solve_galahad_ws, options, inform)
             if (inform%status /= 0) goto 100
          case default
             inform%status = NLLS_ERROR_UNSUPPORTED_METHOD
             GoTo 100
          end select
       end if
    end if

    if (options%scale > 0) then
       call setup_workspace_generate_scaling(n,m,w%generate_scaling_ws,options,inform)
       if (inform%status /= 0) goto 100
    end if
    w%allocated = .true.
100 Continue

    If (inform%status/=0) Call remove_workspace_calculate_step(w,options,tenJ,inner_workspace)

  end subroutine setup_workspace_calculate_step

  recursive subroutine remove_workspace_calculate_step(w,options,tenJ, inner_workspace)
    implicit none
    type( calculate_step_work ), intent(inout) :: w
    type( nlls_options ), intent(in) :: options
    type( tenJ_type), Intent(inout) :: tenJ
    Type( NLLS_workspace), Intent(InOut) :: inner_workspace
    Integer :: ierr_dummy

    if (allocated(w%A)) deallocate(w%A, stat=ierr_dummy)
    if (allocated(w%v)) deallocate(w%v, stat=ierr_dummy)
    if (allocated(w%xxt)) deallocate(w%xxt, stat=ierr_dummy)
    if (allocated(w%extra_scale)) deallocate(w%extra_scale, stat=ierr_dummy)
    if (allocated(w%scale)) deallocate(w%scale, stat=ierr_dummy)

    call remove_workspace_evaluate_model(w%evaluate_model_ws,&
         options)

    if (options%model == 4) then
       call remove_workspace_solve_newton_tensor(&
            w%solve_newton_tensor_ws, &
            options, tenJ, inner_workspace)
    else
       if (w%dogleg_ws%allocated) then
          call remove_workspace_dogleg(w%dogleg_ws, options)
       end if
       if (w%AINT_tr_ws%allocated) then
          call remove_workspace_AINT_tr(w%AINT_tr_ws, options)
       end if
       if (w%more_sorensen_ws%allocated) then
          call remove_workspace_more_sorensen(&
               w%more_sorensen_ws,options)
       end if
       if (w%solve_galahad_ws%allocated) then
          call remove_workspace_solve_galahad(&
               w%solve_galahad_ws, options)
       end if
       if (w%regularization_solver_ws%allocated) then
          call remove_workspace_regularization_solver(&
               w%regularization_solver_ws, options)
       end if
    end if
    if (options%scale > 0) call remove_workspace_generate_scaling(w%generate_scaling_ws,options)

    w%allocated = .false.
  end subroutine remove_workspace_calculate_step

  subroutine setup_workspace_dogleg(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( dogleg_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate(w%d_sd(n),w%d_gn(n),w%ghat(n),w%Jg(m),stat=inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%bad_alloc = 'setup_workspace_dogleg'
      inform%status = NLLS_ERROR_ALLOCATION
      goto 100
    End If
    ! setup space for solve_LLS
    call setup_workspace_solve_LLS(n,m,w%solve_LLS_ws,options,inform)
    if (inform%status /= 0 ) goto 100
    ! setup space for evaluate_model
    call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,options,inform)
    if (inform%status /= 0 ) goto 100

    w%allocated = .true.

100 continue

    If (inform%status /= 0) Call remove_workspace_dogleg(w, options)

  end subroutine setup_workspace_dogleg

  subroutine remove_workspace_dogleg(w,options)
    implicit none
    type( dogleg_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%d_sd )) deallocate(w%d_sd,stat=ierr_dummy )
    if(allocated( w%d_gn )) deallocate(w%d_gn,stat=ierr_dummy )
    if(allocated( w%ghat )) deallocate(w%ghat,stat=ierr_dummy )
    if(allocated( w%Jg )) deallocate(w%Jg,stat=ierr_dummy )

    ! deallocate space for solve_LLS
    call remove_workspace_solve_LLS(w%solve_LLS_ws,options)
    ! deallocate space for evaluate_model
    call remove_workspace_evaluate_model(w%evaluate_model_ws,options)

    w%allocated = .false.
  end subroutine remove_workspace_dogleg

  subroutine setup_workspace_solve_LLS(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( solve_LLS_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform
    integer :: lwork

    inform%status = 0
    lwork = max(1, min(m,n) + max(min(m,n), 1)*4)
    allocate( w%temp(max(m,n)),w%work(lwork),w%Jlls(n*m), stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      Call remove_workspace_solve_LLS(w,options)
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_solve_LLS"
    Else
      w%allocated = .true.
    End If
  end subroutine setup_workspace_solve_LLS

  subroutine remove_workspace_solve_LLS(w,options)
    implicit none
    type( solve_LLS_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%temp )) deallocate( w%temp, stat=ierr_dummy )
    if(allocated( w%work )) deallocate( w%work, stat=ierr_dummy )
    if(allocated( w%Jlls )) deallocate( w%Jlls, stat=ierr_dummy )

    w%allocated = .false.
  end subroutine remove_workspace_solve_LLS

  subroutine setup_workspace_evaluate_model(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( evaluate_model_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate( w%Jd(m),w%dH(n**2),w%dHd(m),w%Hd(n), stat = inform%alloc_status )
    If (inform%alloc_status /= 0) Then
      Call remove_workspace_evaluate_model(w,options)
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = 'evaluate_model'
    Else
      w%allocated = .true.
    End If
  end subroutine setup_workspace_evaluate_model

  subroutine remove_workspace_evaluate_model(w,options)
    implicit none
    type( evaluate_model_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%Jd )) deallocate( w%Jd, stat=ierr_dummy )
    if(allocated( w%dHd )) deallocate( w%dHd, stat=ierr_dummy )
    if(allocated( w%Hd )) deallocate( w%Hd, stat=ierr_dummy )

    w%allocated = .false.
  end subroutine remove_workspace_evaluate_model

  subroutine setup_workspace_AINT_tr(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( AINT_tr_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate(w%B(n,n),w%p0(n),w%p1(n),w%M0(2*n,2*n),w%M1(2*n,2*n),  &
             w%M0_small(n,n),w%M1_small(n,n),w%y(2*n),w%gtg(n,n),   &
             w%q(n),w%LtL(n,n),w%y_hardcase(n,2),w%By_hardcase(n,2),&
             stat = inform%alloc_status)

    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_AINT_tr"
      goto 100
    End If
    ! setup space for max_eig
    call setup_workspace_max_eig(n,m,w%max_eig_ws,options,inform)
    if (inform%status /= 0) goto 100
    call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,options,inform)
    if (inform%status /= 0) goto 100
    ! setup space for the solve routine
    if ((options%model .ne. 1)) then
       call setup_workspace_minus_solve_general(n,m,w%minus_solve_general_ws,options,inform)
       if (inform%status /= 0 ) goto 100
    end if

    w%allocated = .true.
    ! zero-out y_hardcase
    w%y_hardcase(:,:) = 0.0_wp


100 continue

    If (inform%status/=0) Call remove_workspace_AINT_tr(w, options)

  end subroutine setup_workspace_AINT_tr

  subroutine remove_workspace_AINT_tr(w,options)
    implicit none
    type( AINT_tr_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%B )) deallocate(w%B,stat=ierr_dummy)
    if(allocated( w%p0 )) deallocate(w%p0,stat=ierr_dummy)
    if(allocated( w%p1 )) deallocate(w%p1,stat=ierr_dummy)
    if(allocated( w%M0 )) deallocate(w%M0,stat=ierr_dummy)
    if(allocated( w%M1 )) deallocate(w%M1,stat=ierr_dummy)
    if(allocated( w%M0_small )) deallocate(w%M0_small,stat=ierr_dummy)
    if(allocated( w%M1_small )) deallocate(w%M1_small,stat=ierr_dummy)
    if(allocated( w%y )) deallocate(w%y,stat=ierr_dummy)
    if(allocated( w%gtg )) deallocate(w%gtg,stat=ierr_dummy)
    if(allocated( w%q )) deallocate(w%q,stat=ierr_dummy)
    if(allocated( w%LtL )) deallocate(w%LtL,stat=ierr_dummy)
    if(allocated( w%y_hardcase )) deallocate(w%y_hardcase,stat=ierr_dummy)
    if(allocated( w%By_hardcase )) deallocate(w%By_hardcase,stat=ierr_dummy)
    ! setup space for max_eig
    call remove_workspace_max_eig(w%max_eig_ws,options)
    call remove_workspace_evaluate_model(w%evaluate_model_ws,options)
    ! setup space for the solve routine
    if (options%model .ne. 1) then
       call remove_workspace_minus_solve_general(w%minus_solve_general_ws,options)
    end if

    w%allocated = .false.
  end subroutine remove_workspace_AINT_tr

  subroutine setup_workspace_min_eig_symm(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( min_eig_symm_work), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    real(wp), allocatable :: workquery(:)
    integer :: lwork, eigsout, ierr_dummy
    Real(Kind=wp) :: w_dummy(1), z_dummy(1,1)

    inform%status = 0
    allocate(w%A(n,n),workquery(1),stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_min_eig_symm"
      GoTo 100
    End If

    if (options%subproblem_eig_fact) then
       allocate(w%ew(n), stat = inform%alloc_status)
       If (inform%alloc_status /= 0) Then
         inform%status = NLLS_ERROR_ALLOCATION
         inform%bad_alloc = "setup_workspace_min_eig_symm"
         GoTo 100
       End If
       call dsyev('V', & ! both ew's and ev's
            'U', & ! upper triangle of A
            n, w%A, max(1,n), & ! data about A
            w%ew, workquery, -1, &
            inform%external_return)
       If (inform%external_return .ne. 0) then
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = "lapack_dsyev"
          goto 100
       End If
    else
      ! Note: w%ew is allocated with size N to comply with LAPACK interface
       allocate( w%iwork(5*n),w%ifail(n),w%ew(n),stat=inform%alloc_status )
       If (inform%alloc_status /= 0) Then
         inform%status = NLLS_ERROR_ALLOCATION
         inform%bad_alloc = "setup_workspace_min_eig_symm"
         GoTo 100
       End If
       w_dummy(1) = 1.0_wp
       z_dummy(1,1) = 1.0_wp
       ! make a workspace query to dsyevx
       call dsyevx( 'V',& ! get both ew's and ev's
            'I',& ! just the numbered eigenvalues
            'U',& ! upper triangle of A
            n, w%A, n, &
            1.0_wp, 1.0_wp, & ! not used for RANGE = 'I'
            1, 1, & ! only find the first eigenpair
            0.5_wp, & ! abstol for the eigensolver
            eigsout, & ! total number of eigs found
            w_dummy, z_dummy, & ! the eigenvalue and eigenvector
            n, & ! ldz (the eigenvector array)
            workquery, -1, w%iwork, &  ! workspace
            w%ifail, & ! array containing indicies of non-converging ews
            inform%external_return)
       If (inform%external_return .ne. 0) Then
         inform%status = NLLS_ERROR_FROM_EXTERNAL
         inform%external_name = "lapack_dsyevx"
         goto 100
       End If
    end if
    lwork = int(workquery(1))
    if (allocated(workquery)) deallocate(workquery, stat=ierr_dummy)
    allocate( w%work(lwork), stat = inform%alloc_status )
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_min_eig_symm"
      GoTo 100
    End If

    w%allocated = .true.

100 continue

    If (inform%status/=0) Call remove_workspace_min_eig_symm(w,options)

  end subroutine setup_workspace_min_eig_symm

  subroutine remove_workspace_min_eig_symm(w,options)
    implicit none
    type( min_eig_symm_work), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%A )) deallocate(w%A, stat=ierr_dummy)

    if (options%subproblem_eig_fact) then
       if(allocated( w%ew )) deallocate(w%ew, stat=ierr_dummy)
    else
       if(allocated( w%iwork )) deallocate( w%iwork, stat=ierr_dummy )
       if(allocated( w%ifail )) deallocate( w%ifail, stat=ierr_dummy )
       if(allocated( w%ew )) deallocate(w%ew, stat=ierr_dummy)
    end if
    if(allocated( w%work )) deallocate( w%work, stat=ierr_dummy )

    w%allocated = .false.
  end subroutine remove_workspace_min_eig_symm

  subroutine setup_workspace_max_eig(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n, m
    type( max_eig_work), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform), intent(inout) :: inform
    real(wp), allocatable :: workquery(:)
    integer :: lwork, ierr_dummy
    Real(Kind=wp) :: A_dummy(1,1), B_dummy(1,1), vl_dummy(1,1), vr_dummy(1,1)

    inform%status = 0
    allocate( w%alphaR(2*n),w%alphaI(2*n),w%beta(2*n),w%vr(2*n,2*n),           &
      w%ew_array(2*n),workquery(1),w%nullindex(2*n),w%vecisreal(2*n),          &
      stat=inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_max_eig"
      GoTo 100
    End If
    A_dummy(1,1) = 1.0_wp
    B_dummy(1,1) = 1.0_wp
    vl_dummy(1,1) = 0.1_wp
    vr_dummy(1,1) = 0.1_wp
    ! make a workspace query to dggev
    call dggev('N', & ! No left eigenvectors
         'V', &! Yes right eigenvectors
         2*n, A_dummy, max(1,2*n), B_dummy, max(1,2*n), &
         w%alphaR, W%alphaI, w%beta, & ! eigenvalue data
         vl_dummy, max(1,2*n), & ! not referenced
         vr_dummy, max(1,2*n), & ! right eigenvectors
         workquery, -1, inform%external_return)
    If (inform%external_return > 0) Then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = "lapack_dggev"
      GoTo 100
    End If
    lwork = int(workquery(1))
    w%alphaR = 0.0_wp
    w%alphaI = 0.0_wp
    w%beta = 0.0_wp
    w%ew_array = 0.0_wp

    if (allocated(workquery)) deallocate(workquery,stat=ierr_dummy)
    allocate( w%work(lwork), stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_max_eig"
      GoTo 100
    End If

    w%allocated = .true.
100 Continue

    If (inform%status/=0) Call remove_workspace_max_eig(w,options)

  end subroutine setup_workspace_max_eig

  subroutine remove_workspace_max_eig(w,options)
    implicit none
    type( max_eig_work), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%alphaR )) deallocate( w%alphaR, stat=ierr_dummy )
    if(allocated( w%alphaI )) deallocate( w%alphaI, stat=ierr_dummy )
    if(allocated( w%beta )) deallocate( w%beta, stat=ierr_dummy )
    if(allocated( w%vr )) deallocate( w%vr, stat=ierr_dummy )
    if(allocated( w%ew_array )) deallocate( w%ew_array, stat=ierr_dummy )
    if(allocated( w%work )) deallocate( w%work, stat=ierr_dummy )
    if(allocated( w%nullindex )) deallocate( w%nullindex, stat=ierr_dummy )
    if(allocated( w%vecisreal )) deallocate( w%vecisreal, stat=ierr_dummy )

    w%allocated = .false.
  end subroutine remove_workspace_max_eig

  subroutine setup_workspace_minus_solve_general(n, m, w, options, inform)
    implicit none
    integer, intent(in) :: n, m
    type( minus_solve_general_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform), intent(inout) :: inform

    inform%status = 0
    allocate( w%A(n,n),w%ipiv(n),stat=inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      Call remove_workspace_minus_solve_general(w,options)
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_minus_solve_general"
    Else
      w%allocated = .true.
    End If
  end subroutine setup_workspace_minus_solve_general

  subroutine remove_workspace_minus_solve_general(w, options)
    implicit none
    type( minus_solve_general_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%A )) deallocate( w%A, stat=ierr_dummy )
    if(allocated( w%ipiv )) deallocate( w%ipiv, stat=ierr_dummy )

    w%allocated = .false.
  end subroutine remove_workspace_minus_solve_general

  subroutine setup_workspace_solve_galahad(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n,m
    type( solve_galahad_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate(w%ev(n,n),w%v_trans(n),w%ew(n),w%d_trans(n),w%scale_c(n),         &
      w%scale_h(n),stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_solve_galahad"
      GoTo 100
    End If

    call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,options,inform)
    if (inform%status /= 0) goto 100

    w%allocated = .true.

100 continue

    ! This call duplicates some deallocs but it is ok
    If (inform%status/=0) Call remove_workspace_solve_galahad(w,options)

  end subroutine setup_workspace_solve_galahad

  subroutine remove_workspace_solve_galahad(w,options)
    implicit none
    type( solve_galahad_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%ev )) deallocate(w%ev, stat=ierr_dummy)
    if(allocated( w%v_trans )) deallocate(w%v_trans, stat=ierr_dummy)
    if(allocated( w%ew )) deallocate(w%ew, stat=ierr_dummy)
    if(allocated( w%d_trans )) deallocate(w%d_trans, stat=ierr_dummy)
    if(allocated( w%scale_c )) deallocate(w%scale_c, stat=ierr_dummy)
    if(allocated( w%scale_h )) deallocate(w%scale_h, stat=ierr_dummy)

    call remove_workspace_all_eig_symm(w%all_eig_symm_ws,options)

    w%allocated = .false.
  end subroutine remove_workspace_solve_galahad

  subroutine setup_workspace_regularization_solver(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n,m
    type ( regularization_solver_work ), INTENT( INOUT) :: w
    type ( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate(w%LtL(n,n),w%AplusSigma(n,n),stat = inform%alloc_status)
    If (inform%alloc_status /= 0) Then
      Call remove_workspace_regularization_solver(w,options)
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_regularization_solver"
    Else
      w%allocated = .true.
    End If
  end subroutine setup_workspace_regularization_solver

  subroutine remove_workspace_regularization_solver(w,options)
    implicit none
    type( regularization_solver_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%LtL )) deallocate(w%LtL, stat=ierr_dummy)
    if(allocated( w%AplusSigma )) deallocate(w%AplusSigma, stat=ierr_dummy)

    w%allocated = .false.
  end subroutine remove_workspace_regularization_solver

  subroutine setup_workspace_all_eig_symm(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n,m
    type( all_eig_symm_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    real(wp), allocatable :: workquery(:)
    real(wp) :: A_dummy(1)
    integer :: lwork, ierr_dummy

    inform%status = 0
    If (allocated(w%ew)) Deallocate(w%ew, stat=ierr_dummy)
    allocate(workquery(1),w%ew(n), stat = inform%alloc_status)
    If (inform%alloc_status /= 0 ) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_all_eig_symm"
      GoTo 100
    End If
    A_dummy = 1.0_wp
    w%ew(:) = 1.0_wp
    call dsyev('V', & ! both ew's and ev's
         'U', & ! upper triangle of A
         n, A_dummy, max(1,n), & ! data about A
         w%ew, workquery, -1, &
         inform%external_return)
    If (inform%external_return .ne. 0) Then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = "lapack_dsyev"
      GoTo 100
    End If

    lwork = int(workquery(1))
    if (allocated(workquery)) deallocate(workquery, stat=ierr_dummy)
    allocate( w%work(lwork), stat = inform%alloc_status )
    If (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_all_eig_symm"
      GoTo 100
    End If

    w%allocated = .true.

100 continue
    If (allocated(w%ew)) Deallocate(w%ew, stat=ierr_dummy)
    If (inform%status/=0) Call remove_workspace_all_eig_symm(w,options)

  end subroutine setup_workspace_all_eig_symm

  subroutine remove_workspace_all_eig_symm(w,options)
    implicit none
    type( all_eig_symm_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%work )) deallocate( w%work, stat=ierr_dummy )
    w%allocated = .false.
  end subroutine remove_workspace_all_eig_symm

  subroutine setup_workspace_more_sorensen(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n,m
    type( more_sorensen_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    allocate(w%LtL(n,n),w%q(n),w%y1(n),w%AplusSigma(n,n),w%norm_work(n),       &
      stat = inform%alloc_status)
    if (inform%alloc_status /= 0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_more_sorensen"
      goto 100
    End If

    call setup_workspace_min_eig_symm(n,m,w%min_eig_symm_ws,options,inform)
    if (inform%status /= 0) goto 100

    w%allocated = .true.

100 continue

    If (inform%status/=0) Call remove_workspace_more_sorensen(w,options)

  end subroutine setup_workspace_more_sorensen

  subroutine remove_workspace_more_sorensen(w,options)
    implicit none
    type( more_sorensen_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%LtL )) deallocate(w%LtL, stat=ierr_dummy)
    if(allocated( w%q )) deallocate(w%q, stat=ierr_dummy)
    if(allocated( w%y1 )) deallocate(w%y1, stat=ierr_dummy)
    if(allocated( w%AplusSigma )) deallocate(w%AplusSigma, stat=ierr_dummy)
    if(allocated( w%norm_work )) deallocate(w%norm_work, stat=ierr_dummy)

    call remove_workspace_min_eig_symm(w%min_eig_symm_ws,options)
    w%allocated = .false.
  end subroutine remove_workspace_more_sorensen


  subroutine setup_workspace_generate_scaling(n,m,w,options,inform)
    implicit none
    integer, intent(in) :: n,m
    type( generate_scaling_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(inout) :: inform

    inform%status = 0
    If (options%scale==4) Then
      allocate(w%tempvec(n),w%diag(n),w%ev(n,n),stat=inform%alloc_status)
    Else
      allocate(w%diag(n),w%ev(n,n),stat=inform%alloc_status)
    End If
    If (inform%alloc_status/=0) Then
      inform%status = NLLS_ERROR_ALLOCATION
      inform%bad_alloc = "setup_workspace_generate_scaling"
      GoTo 100
    End If

    if (options%scale==4) then
       call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,options,inform)
       if (inform%status/=0) goto 100
    end if

    w%diag = 1.0_wp
    w%allocated = .true.

100 continue

    If (inform%status/=0) Call remove_workspace_generate_scaling(w,options)

  end subroutine setup_workspace_generate_scaling

  subroutine remove_workspace_generate_scaling(w,options)
    implicit none
    type( generate_scaling_work ), INTENT( INOUT) :: w
    type( nlls_options ), intent(in) :: options
    Integer :: ierr_dummy

    if(allocated( w%diag )) deallocate( w%diag, stat=ierr_dummy )
    if(allocated( w%ev )) deallocate( w%ev, stat=ierr_dummy )

    if (options%scale==4) then
       call remove_workspace_all_eig_symm(w%all_eig_symm_ws,options)
       if(allocated( w%tempvec )) deallocate( w%tempvec, stat=ierr_dummy )
    end if

    w%allocated = .false.
  end subroutine remove_workspace_generate_scaling

  subroutine setup_bounds_type(w, n, lower_bounds, upper_bounds, options, inform)
    Implicit None
    type(box_type), Intent(InOut)             :: w
    Integer, Intent(In)                       :: n
    Real(Kind=wp), Intent(In), Optional       :: lower_bounds(n), upper_bounds(n)
    Type(NLLS_options), Intent(In)            :: options
    Type(NLLS_inform), Intent(InOut)          :: inform
    
    Integer                                   :: i, ierr
    Logical                                   :: has_box, lower_usr, upper_usr
    Real(wp)                                  :: lower, upper

    has_box = .False.
    w%has_box = has_box
    w%prjchd = .False.
    w%quad_i = 0
    w%quad_c = 0.0_wp
    w%quad_q = 0.0_wp
    lower_usr = Present(lower_bounds)
    upper_usr = Present(upper_bounds)

    if ( (.Not. lower_usr) .And. (.Not. upper_usr) ) Then
      Go To 100
    end if
    
    Do i = 1, n
        If (lower_usr) Then
          lower = lower_bounds(i)
        Else
          lower = -options%box_bigbnd
        End if
        If (upper_usr) Then
          upper = upper_bounds(i)
        Else
          upper = options%box_bigbnd
        End if

       If ( lower <= upper .And. lower == lower .And. upper == upper ) Then
          If (-options%box_bigbnd < lower .Or. upper < options%box_bigbnd) Then
             has_box = .True.
          End If
       Else
          inform%status = NLLS_ERROR_BAD_BOX_BOUNDS
          Go To 100
       End If
    End Do

    w%has_box = has_box

    If (has_box) Then
      ! Clear all arrays...
      If (allocated(w%blx)) deallocate(w%blx, Stat=ierr)
      If (allocated(w%bux)) deallocate(w%bux, Stat=ierr)
      If (allocated(w%pdir)) deallocate(w%pdir, Stat=ierr)
      If (allocated(w%normFref)) deallocate(w%normFref, Stat=ierr)
      If (allocated(w%sk)) deallocate(w%sk, Stat=ierr)
      If (allocated(w%g)) deallocate(w%g, Stat=ierr)
      Allocate(w%blx(n), w%bux(n), w%pdir(n), w%g(n),      &
        w%normFref(options%box_nFref_max), w%sk(n), Stat=ierr)
      if (ierr /= 0) Then
        If (allocated(w%blx)) deallocate(w%blx, Stat=ierr)
        If (allocated(w%bux)) deallocate(w%bux, Stat=ierr)
        If (allocated(w%pdir)) deallocate(w%pdir, Stat=ierr)
        If (allocated(w%normFref)) deallocate(w%normFref, Stat=ierr)
        If (allocated(w%sk)) deallocate(w%sk, Stat=ierr)
        If (allocated(w%g)) deallocate(w%g, Stat=ierr)
        inform%status = NLLS_ERROR_ALLOCATION
        inform%bad_alloc = 'ral_nlls_box'
        Go To 100
      end if
      w%normfref(1:options%box_nFref_max) = -1.0e-20_wp
      Do i = 1, n
        If (lower_usr) Then
          lower = lower_bounds(i)
        Else
          lower = -options%box_bigbnd
        End if
        If (upper_usr) Then
          upper = upper_bounds(i)
        Else
          upper = options%box_bigbnd
        End if
        w%blx(i) = max(lower, -options%box_bigbnd)
        w%bux(i) = min(upper, options%box_bigbnd)
      End Do
    End If

100 Continue

  end subroutine setup_bounds_type

  Subroutine remove_workspace_bounds(w)

    Implicit None
    type(box_type), Intent(InOut)             :: w
    Integer                                   :: ierr_dummy
    Continue

    If (allocated(w%blx)) deallocate(w%blx, Stat=ierr_dummy)
    If (allocated(w%bux)) deallocate(w%bux, Stat=ierr_dummy)
    If (allocated(w%pdir)) deallocate(w%pdir, Stat=ierr_dummy)
    If (allocated(w%normFref)) deallocate(w%normFref, Stat=ierr_dummy)
    If (allocated(w%sk)) deallocate(w%sk, Stat=ierr_dummy)
    If (allocated(w%g)) deallocate(w%g, Stat=ierr_dummy)

  End Subroutine remove_workspace_bounds



end module ral_nlls_workspaces
