! ral_nlls_workspaces :: module to keep all the workspaces

module ral_nlls_workspaces 

  implicit none 
  
  private
  
  ! define derived types and subroutines for workspace arrays.
  
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
  real (kind = wp), parameter :: three = 3.0
  real (kind = wp), parameter :: half = 0.5
  real (kind = wp), parameter :: sixteenth = 0.0625

  
  
    TYPE, PUBLIC :: NLLS_options
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

!  which linear least squares solver should we use?
     
     INTEGER :: lls_solver = 1
        
!   overall convergence tolerances. The iteration will terminate when the
!   norm of the gradient of the objective function is smaller than 
!      MAX( %stop_g_absolute, %stop_g_relative * norm of the initial gradient)
!   or if the norm of objective function is smaller than 
!      MAX( %stop%f_absolute, %stop_f_relative * initial norm of the function)
!     or if the step is less than %stop_s

     REAL ( KIND = wp ) :: stop_g_absolute = tenm5
     REAL ( KIND = wp ) :: stop_g_relative = tenm8
     REAL ( KIND = wp ) :: stop_f_absolute = tenm8
     REAL ( KIND = wp ) :: stop_f_relative = tenm8
     REAL ( KIND = wp ) :: stop_s = epsmch

     
!   should we scale the initial trust region radius?
     
     integer :: relative_tr_radius = 0

!   if relative_tr_radius == 1, then pick a scaling parameter
!   Madsen, Nielsen and Tingleff say pick this to be 1e-6, say, if x_0 is good,
!   otherwise 1e-3 or even 1 would be good starts...
     
     real (kind = wp) :: initial_radius_scale = 1.0

!   if relative_tr_radius /= 1, then set the 
!   initial value for the trust-region radius (-ve => ||g_0||)
     
     REAL ( KIND = wp ) :: initial_radius = hundred
     
!   for the newton tensor model, allow a base tr raidius to allow an inherent
!   regularization in the problem that can't be changed
!   ( so we minimize 0.5 * (\sum (f_i)^2 + sigma_k ||s||^2) ), using another reg parameter
!   on top of this
!   (undocumented control variable)
     
     REAL ( KIND = wp ) :: base_regularization = zero

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
     REAL ( KIND = wp ) :: regularization_term = 1e-2_wp
     REAL ( KIND = wp ) :: regularization_power = 2.0_wp
     REAL ( KIND = wp ) :: regularization_weight = 1.0_wp

     

!   maximum permitted trust-region radius

     REAL ( KIND = wp ) :: maximum_radius = ten ** 8

!   a potential iterate will only be accepted if the actual decrease
!    f - f(x_new) is larger than %eta_successful times that predited
!    by a quadratic model of the decrease. The trust-region radius will be
!    increased if this relative decrease is greater than %eta_very_successful
!    but smaller than %eta_too_successful

     REAL ( KIND = wp ) :: eta_successful = ten ** ( - 8 )! ten ** ( - 8 ) 
     REAL ( KIND = wp ) :: eta_success_but_reduce = ten ** ( - 8 ) !0.25_wp
     REAL ( KIND = wp ) :: eta_very_successful = point9!0.75_wp!point9 
     REAL ( KIND = wp ) :: eta_too_successful = two

!   on very successful iterations, the trust-region radius will be increased by
!    the factor %radius_increase, while if the iteration is unsuccessful, the 
!    radius will be decreased by a factor %radius_reduce but no more than
!    %radius_reduce_max

     REAL ( KIND = wp ) :: radius_increase = two
     REAL ( KIND = wp ) :: radius_reduce = half
     REAL ( KIND = wp ) :: radius_reduce_max = sixteenth

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
     

!   scale the variables?
!   0 - no scaling
!   1 - use the scaling in GSL (W s.t. W_ii = ||J(i,:)||_2^2)
!       tiny values get set to one       
!   2 - scale using the approx to the Hessian (W s.t. W = ||H(i,:)||_2^2
     INTEGER :: scale = 1
     REAL(wp) :: scale_max = 1e11
     REAL(wp) :: scale_min = 1e-11
     LOGICAL :: scale_trim_min = .TRUE.
     LOGICAL :: scale_trim_max = .TRUE.
     LOGICAL :: scale_require_increase = .FALSE.

     logical :: calculate_svd_J = .false.

     logical :: setup_workspaces = .true.
     logical :: remove_workspaces = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! M O R E - S O R E N S E N   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

    integer  :: more_sorensen_maxits = 500
    real(wp) :: more_sorensen_shift = 1e-13
    real(wp) :: more_sorensen_tiny = 10.0_wp * epsmch
    real(wp) :: more_sorensen_tol = 1e-3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! H Y B R I D   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what's the tolerance such that ||J^T f || < tol * 0.5 ||f||^2 triggers a switch
    real(wp) :: hybrid_tol = 2.0!0.02

! how many successive iterations does the above condition need to hold before we switch?
    integer  :: hybrid_switch_its = 1!3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! T E N S O R   M O D E L   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what regularization should we use?
    real(wp) :: reg_order = -one

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
     
     INTEGER :: iter

!  the number of inner iterations performed
     
     INTEGER :: inner_iter = 0
       
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

!  test on the size of the step satisfied?
     
     integer :: convergence_norms = 0
     
!  vector of residuals 
     
     real(wp), allocatable :: resvec(:)

!  vector of gradients 
     
     real(wp), allocatable :: gradvec(:)

!  vector of smallest singular values

     real(wp), allocatable :: smallest_sv(:)
     
!  vector of largest singular values

     real(wp), allocatable :: largest_sv(:)

!  the value of the objective function at the best estimate of the solution 
!   determined by NLLS_solve

     REAL ( KIND = wp ) :: obj = HUGE( one )

!  the norm of the gradient of the objective function at the best estimate 
!   of the solution determined by NLLS_solve

     REAL ( KIND = wp ) :: norm_g = HUGE( one )

! the norm of the gradient, scaled by the norm of the residual
     
     REAL ( KIND = wp ) :: scaled_g = HUGE( one ) 

! error returns from external subroutines 
     
     INTEGER :: external_return = 0 
     
! name of external program that threw and error
     
     CHARACTER ( LEN = 80 ) :: external_name = REPEAT( ' ', 80 )

  END TYPE nlls_inform


  TYPE, public :: NLLS_ERROR
     INTEGER :: MAXITS = -1
     INTEGER :: EVALUATION = -2
     INTEGER :: UNSUPPORTED_MODEL = -3
     INTEGER :: FROM_EXTERNAL = -4
     INTEGER :: UNSUPPORTED_METHOD = -5
     INTEGER :: ALLOCATION = -6
     INTEGER :: MAX_TR_REDUCTIONS = -7
     INTEGER :: X_NO_PROGRESS = -8
     INTEGER :: N_GT_M = -9  
     INTEGER :: BAD_TR_STRATEGY = -10 
     INTEGER :: FIND_BETA = -11 
     INTEGER :: BAD_SCALING = -12 
     INTEGER :: WORKSPACE_ERROR = -13
     ! dogleg errors
     INTEGER :: DOGLEG_MODEL = -101
     ! AINT errors
     INTEGER :: AINT_EIG_IMAG = -201
     INTEGER :: AINT_EIG_ODD = -202
     ! More-Sorensen errors
     INTEGER :: MS_MAXITS = -301
     INTEGER :: MS_TOO_MANY_SHIFTS = -302
     INTEGER :: MS_NO_PROGRESS = -303
     ! DTRS errors
     ! Tensor model errors
     INTEGER :: NO_SECOND_DERIVATIVES = -401
     INTEGER :: NT_BAD_SUBPROBLEM = -402

  END TYPE NLLS_ERROR

  type(NLLS_ERROR), public :: ERROR


  type, public :: params_base_type
     ! deliberately empty
  end type params_base_type


  type, extends( params_base_type ), public :: tensor_params_type
     ! blank?
     real(wp), dimension(:), allocatable :: f
     real(wp), dimension(:), allocatable :: J
     real(wp), dimension(:,:,:), allocatable :: Hi
     real(wp) :: Delta
     real(wp) :: p
     integer :: m
     integer :: m1 = 0
     integer :: extra = 0
  end type tensor_params_type


  type, public :: max_eig_work ! workspace for subroutine max_eig
     logical :: allocated = .false.
     real(wp), allocatable :: alphaR(:), alphaI(:), beta(:), vr(:,:)
     real(wp), allocatable :: work(:), ew_array(:)
     integer, allocatable :: nullindex(:)
     logical, allocatable :: vecisreal(:)
     integer :: nullevs_cols
     real(wp), allocatable :: nullevs(:,:)
  end type max_eig_work

  type, public :: solve_general_work ! workspace for subroutine solve_general
     logical :: allocated = .false.
     real(wp), allocatable :: A(:,:)
     integer, allocatable :: ipiv(:)
  end type solve_general_work

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
  end type all_eig_symm_work

  type, public :: apply_scaling_work ! workspace for subrouine apply_scaling
     logical :: allocated = .false.
     real(wp), allocatable :: diag(:)
     real(wp), allocatable :: ev(:,:)
     real(wp), allocatable :: tempvec(:)
     type( all_eig_symm_work ) :: all_eig_symm_ws
  end type apply_scaling_work

  type, public :: solve_newton_tensor_work ! a workspace for solve_newton_tensor
     logical :: allocated = .false.
     real(wp), allocatable :: model_tensor(:)
     type( tensor_params_type ) :: tparams
     type( nlls_options ) :: tensor_options
     integer :: m_in
  end type solve_newton_tensor_work

  type, public :: solve_galahad_work ! workspace for subroutine dtrs_work
     logical :: allocated = .false.
     real(wp), allocatable ::ev(:,:), ew(:), v_trans(:), d_trans(:)
     real(wp) :: reg_order
     type( all_eig_symm_work ) :: all_eig_symm_ws
  end type solve_galahad_work

  type, public :: more_sorensen_work ! workspace for subroutine more_sorensen
     logical :: allocated = .false.
     real(wp), allocatable :: LtL(:,:), AplusSigma(:,:)
     real(wp), allocatable :: q(:), y1(:)
     type( min_eig_symm_work ) :: min_eig_symm_ws
  end type more_sorensen_work

  type, public :: AINT_tr_work ! workspace for subroutine AINT_tr
     logical :: allocated = .false.
     type( max_eig_work ) :: max_eig_ws
     type( evaluate_model_work ) :: evaluate_model_ws
     type( solve_general_work ) :: solve_general_ws
     !       type( solve_spd_work ) :: solve_spd_ws
     REAL(wp), allocatable :: LtL(:,:), B(:,:), p0(:), p1(:)
     REAL(wp), allocatable :: M0(:,:), M1(:,:), y(:), gtg(:,:), q(:)
     REAL(wp), allocatable :: M0_small(:,:), M1_small(:,:)
     REAL(wp), allocatable :: y_hardcase(:,:)
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
     real(wp), allocatable :: v(:), extra_scale(:)
     type( AINT_tr_work ) :: AINT_tr_ws
     type( dogleg_work ) :: dogleg_ws
     type( solve_newton_tensor_work ) :: solve_newton_tensor_ws
     type( more_sorensen_work ) :: more_sorensen_ws
     type( solve_galahad_work ) :: solve_galahad_ws
     type( evaluate_model_work) :: evaluate_model_ws
     type( apply_scaling_work ) :: apply_scaling_ws
  end type calculate_step_work

  type, public :: get_svd_J_work ! workspace for subroutine get_svd_J
     logical :: allocated = .false.
     real(wp), allocatable :: Jcopy(:), S(:), work(:)
  end type get_svd_J_work

  type, public :: NLLS_workspace ! all workspaces called from the top level
     logical :: allocated = .false.
     integer :: first_call = 1
     integer :: iter = 0 
     real(wp) :: normF0, normJF0, normF, normJF
     real(wp) :: normJFold, normJF_Newton
     real(wp) :: Delta
     real(wp) :: normd
     real(wp) :: reg_order = two ! reg. by + 1/p || \sigma || ** p
     logical :: use_second_derivatives = .false.
     integer :: hybrid_count = 0
     real(wp) :: hybrid_tol = 1.0
     real(wp), allocatable :: fNewton(:), JNewton(:), XNewton(:)
     real(wp), allocatable :: J(:)
     real(wp), allocatable :: f(:), fnew(:)
     real(wp), allocatable :: hf(:), hf_temp(:)
     real(wp), allocatable :: d(:), g(:), Xnew(:)
     real(wp), allocatable :: y(:), y_sharp(:), g_old(:), g_mixed(:)
     real(wp), allocatable :: ysharpSks(:), Sks(:)
     real(wp), allocatable :: resvec(:), gradvec(:)
     real(wp), allocatable :: largest_sv(:), smallest_sv(:)
     type ( calculate_step_work ) :: calculate_step_ws
     type ( get_svd_J_work ) :: get_svd_J_ws
     real(wp) :: tr_nu = 2.0
     integer :: tr_p = 3
  end type NLLS_workspace

  type, public :: tenJ_type ! workspace for evaltensor_J
     logical :: allocated = .false.
     real(wp), allocatable :: Hs(:), Js(:) ! work arrays for evaltensor_f
  end type tenJ_type
  type( tenJ_type ), public :: tenJ
  type( nlls_workspace ), public :: inner_workspace ! to be used to solve recursively    

  public :: setup_workspaces, remove_workspaces

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                       !!
  !! W O R K S P A C E   S E T U P   S U B R O U T I N E S !!
  !!                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine setup_workspaces(workspace,n,m,options,inform)

    type( NLLS_workspace ), intent(inout) :: workspace
    type( nlls_options ), intent(in) :: options
    integer, intent(in) :: n,m
    type( nlls_inform ), intent(out) :: inform

    if (.not. allocated(workspace%y)) then
       allocate(workspace%y(n), stat = inform%alloc_status)
       if (inform%alloc_status .ne. 0) goto 1000
       workspace%y = zero
    end if
    if (.not. allocated(workspace%y_sharp)) then
       allocate(workspace%y_sharp(n), stat = inform%alloc_status)
       if (inform%alloc_status .ne. 0) goto 1000
       workspace%y_sharp = zero
    end if

    if (.not. options%exact_second_derivatives) then
       if (.not. allocated(workspace%g_old)) then
          allocate(workspace%g_old(n), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       if (.not. allocated(workspace%g_mixed)) then
          allocate(workspace%g_mixed(n), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       if (.not. allocated(workspace%Sks)) then
          allocate(workspace%Sks(n), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       if (.not. allocated(workspace%ysharpSks)) then
          allocate(workspace%ysharpSks(n), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if

    end if

    if( options%output_progress_vectors ) then 
       if (.not. allocated(workspace%resvec)) then
          allocate(workspace%resvec(options%maxit+1), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       if (.not. allocated(workspace%gradvec)) then
          allocate(workspace%gradvec(options%maxit+1), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
    end if

    if( options%calculate_svd_J ) then
       if (.not. allocated(workspace%largest_sv)) then
          allocate(workspace%largest_sv(options%maxit + 1), &
               stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       if (.not. allocated(workspace%smallest_sv)) then
          allocate(workspace%smallest_sv(options%maxit + 1), &
               stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
       call setup_workspace_get_svd_J(n,m,workspace%get_svd_J_ws, & 
            options, inform)
       if (inform%alloc_status > 0) goto 1010
    end if

    if( .not. allocated(workspace%J)) then
       allocate(workspace%J(n*m), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( .not. allocated(workspace%f)) then
       allocate(workspace%f(m), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( .not. allocated(workspace%fnew)) then 
       allocate(workspace%fnew(m), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( .not. allocated(workspace%hf)) then
       allocate(workspace%hf(n*n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( options%model == 3 ) then
       if( .not. allocated(workspace%hf_temp)) then 
          allocate(workspace%hf_temp(n*n), stat = inform%alloc_status)
          if (inform%alloc_status > 0) goto 1000
       end if
    end if
    if( .not. allocated(workspace%d)) then
       allocate(workspace%d(n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( .not. allocated(workspace%g)) then
       allocate(workspace%g(n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if
    if( .not. allocated(workspace%Xnew)) then
       allocate(workspace%Xnew(n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 1000
    end if

    call setup_workspace_calculate_step(n,m,workspace%calculate_step_ws, & 
         options, inform)
    if (inform%alloc_status > 0) goto 1010

    workspace%allocated = .true.

    return

    ! Error statements
1000 continue ! bad allocation from this subroutine
    inform%bad_alloc = 'setup_workspaces'
    inform%status = ERROR%ALLOCATION
    return

1010 continue ! bad allocation from called subroutine
    return

  end subroutine setup_workspaces

  recursive subroutine remove_workspaces(workspace,options)

    type( NLLS_workspace ), intent(inout) :: workspace
    type( nlls_options ), intent(in) :: options

    if(allocated(workspace%y)) deallocate(workspace%y)
    if(allocated(workspace%y_sharp)) deallocate(workspace%y_sharp)
    if(allocated(workspace%g_old)) deallocate(workspace%g_old)
    if(allocated(workspace%g_mixed)) deallocate(workspace%g_mixed)
    if(allocated(workspace%Sks)) deallocate(workspace%Sks)
    if(allocated(workspace%ysharpSks)) deallocate(workspace%ysharpSks)

    if(allocated(workspace%resvec)) deallocate(workspace%resvec)
    if(allocated(workspace%gradvec)) deallocate(workspace%gradvec)

    if(allocated(workspace%largest_sv)) deallocate(workspace%largest_sv)
    if(allocated(workspace%smallest_sv)) deallocate(workspace%smallest_sv)

    if(allocated(workspace%fNewton)) deallocate(workspace%fNewton )
    if(allocated(workspace%JNewton)) deallocate(workspace%JNewton )
    if(allocated(workspace%XNewton)) deallocate(workspace%XNewton ) 

    if( options%calculate_svd_J ) then
       if (allocated(workspace%largest_sv)) deallocate(workspace%largest_sv)
       if (allocated(workspace%smallest_sv)) deallocate(workspace%smallest_sv)
       call remove_workspace_get_svd_J(workspace%get_svd_J_ws, options)
    end if

    if(allocated(workspace%J)) deallocate(workspace%J ) 
    if(allocated(workspace%f)) deallocate(workspace%f ) 
    if(allocated(workspace%fnew)) deallocate(workspace%fnew ) 
    if(allocated(workspace%hf)) deallocate(workspace%hf ) 
    if(allocated(workspace%hf_temp)) deallocate(workspace%hf_temp) 
    if(allocated(workspace%d)) deallocate(workspace%d ) 
    if(allocated(workspace%g)) deallocate(workspace%g ) 
    if(allocated(workspace%Xnew)) deallocate(workspace%Xnew ) 

    call remove_workspace_calculate_step(workspace%calculate_step_ws,&
         options)

    !       ! evaluate model in the main routine...       
    !       call remove_workspace_evaluate_model(workspace%evaluate_model_ws,options)

    workspace%allocated = .false.

    return

  end subroutine remove_workspaces

  subroutine setup_workspace_get_svd_J(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( get_svd_J_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    integer :: lwork
    character :: jobu(1), jobvt(1)

    allocate( w%Jcopy(n*m), stat = inform%alloc_status)
    if (inform%alloc_status .ne. 0) goto 9000

    allocate( w%S(n), stat = inform%alloc_status)
    if (inform%alloc_status .ne. 0) goto 9000

    allocate(w%work(1))
    jobu  = 'N' ! calculate no left singular vectors
    jobvt = 'N' ! calculate no right singular vectors

    ! make a workspace query....
    call dgesvd( jobu, jobvt, n, m, w%Jcopy, n, w%S, w%S, 1, w%S, 1, & 
         w%work, -1, inform%alloc_status )
    if ( inform%alloc_status .ne. 0 ) goto 9000

    lwork = int(w%work(1))
    if (allocated(w%work)) deallocate(w%work)
    allocate(w%work(lwork))     

    w%allocated = .true.

    return

    ! Error statements
9000 continue  ! bad allocations in this subroutine
    inform%bad_alloc = 'setup_workspace_get_svd_J'
    inform%status = ERROR%ALLOCATION
    return

  end subroutine setup_workspace_get_svd_J

  subroutine remove_workspace_get_svd_J(w,options)
    type( get_svd_J_work ) :: w
    type( nlls_options ) :: options

    if( allocated(w%Jcopy) ) deallocate( w%Jcopy )
    if( allocated(w%S) ) deallocate( w%S )
    if( allocated(w%work) ) deallocate( w%work )

    w%allocated = .false.

    return

  end subroutine remove_workspace_get_svd_J

  subroutine setup_workspace_solve_newton_tensor(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( solve_newton_tensor_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%model_tensor(m), stat=inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(w%tparams%f(m), stat=inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(w%tparams%J(n*m), stat=inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(w%tparams%Hi(n,n,m), stat=inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    w%tparams%m = m

    ! copy options from those input
    w%tensor_options = options
    ! use a hybrid method for the inner loop
    w%tensor_options%model = 2
    w%tensor_options%maxit = 5000
    w%tensor_options%reg_order = -one
    w%tensor_options%output_progress_vectors = .false.

    allocate(tenJ%Hs(n), tenJ%Js(m), stat=inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

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
    case (2)
       w%tensor_options%model = 1
       w%tensor_options%type_of_method = 1 ! make this changable by the user
       w%tensor_options%nlls_method = 4 
       !         w%tensor_options%radius_increase = 2.0_wp
       !         w%tensor_options%radius_reduce = 0.5_wp
       !          w%tensor_options%stop_g_absolute = 1e-14
       !          w%tensor_options%stop_g_relative = 1e-14
       w%tparams%m1 = m
       w%m_in = m + n
    case (3,4,5,6,7)
       w%tensor_options%model = 1
       w%tensor_options%type_of_method = 1
       w%tensor_options%nlls_method = 4
       !         w%tensor_options%stop_g_absolute = 1e-14
       !         w%tensor_options%stop_g_relative = 1e-14
       !          w%tensor_options%radius_increase = 2.0_wp
       !          w%tensor_options%radius_reduce = 0.5_wp
       w%tparams%m1 = m
       w%m_in = m 
       if (options%inner_method == 3 .or. options%inner_method == 5) then 
          w%tensor_options%regularization = 1
          w%tensor_options%regularization_power = 2.0_wp
       else
          w%tensor_options%regularization = 2
          w%tensor_options%regularization_power = 3.0_wp
       end if
       if (options%inner_method == 5 .or. options%inner_method == 6) then 
          w%tensor_options%update_lower_order = .false.
       end if


    end select

    ! setup/remove workspaces manually....
    w%tensor_options%remove_workspaces = .false.
    w%tensor_options%setup_workspaces = .false.
    call setup_workspaces(inner_workspace, n, w%m_in, w%tensor_options, inform)
    if (inform%alloc_status > 0) goto 9000       

    w%allocated = .true.

    return

    ! Error statements
9000 continue  ! bad allocations in this subroutine
    inform%bad_alloc = 'setup_workspace_dogleg'
    inform%status = ERROR%ALLOCATION
    return

  end subroutine setup_workspace_solve_newton_tensor

  subroutine remove_workspace_solve_newton_tensor(w,options)
    type( solve_newton_tensor_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated(w%model_tensor)) deallocate(w%model_tensor)
    if(allocated(w%tparams%f)) deallocate(w%tparams%f)
    if(allocated(w%tparams%J)) deallocate(w%tparams%J)
    if(allocated(w%tparams%Hi)) deallocate(w%tparams%Hi)
    if(allocated(tenJ%Hs)) deallocate(tenJ%Hs)
    if(allocated(tenJ%Js)) deallocate(tenJ%Js)

    call remove_workspaces(inner_workspace, w%tensor_options)

    w%allocated = .false.

    return

  end subroutine remove_workspace_solve_newton_tensor

  recursive subroutine setup_workspace_calculate_step(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( calculate_step_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%A(n,n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(w%v(n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(w%extra_scale(n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000      


    call setup_workspace_evaluate_model(n,m,& 
         w%evaluate_model_ws,options,inform)
    if (inform%alloc_status > 0) goto 9010

    if ( options%model == 4 ) then

       call setup_workspace_solve_newton_tensor(n,m,&
            w%solve_newton_tensor_ws,&
            options, inform)

    else

       select case (options%nlls_method)

       case (1) ! use the dogleg method
          call setup_workspace_dogleg(n,m,w%dogleg_ws, & 
               options, inform)
          if (inform%alloc_status > 0) goto 9010

       case(2) ! use the AINT method
          call setup_workspace_AINT_tr(n,m,w%AINT_tr_ws, & 
               options, inform)
          if (inform%alloc_status > 0) goto 9010

       case(3) ! More-Sorensen 
          call setup_workspace_more_sorensen(n,m,&
               w%more_sorensen_ws,options,inform)
          if (inform%alloc_status > 0) goto 9010

       case (4) ! dtrs (Galahad)
          call setup_workspace_solve_galahad(n,m, & 
               w%solve_galahad_ws, options, inform)
          if (inform%alloc_status > 0) goto 9010
       end select
    end if

    if (options%scale > 0) then
       call setup_workspace_apply_scaling(n,m,w%apply_scaling_ws,options,inform)
       if (inform%status > 0) goto 9010
    end if



    w%allocated = .true.

    return

    ! Error statements
9000 continue  ! bad allocations in this subroutine
    inform%bad_alloc = 'setup_workspace_calculate_step'
    inform%status = ERROR%ALLOCATION
    return

9010 continue  ! bad allocations from dependent subroutine
    return
  end subroutine setup_workspace_calculate_step

  recursive subroutine remove_workspace_calculate_step(w,options)
    type( calculate_step_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options

    if (allocated(w%A)) deallocate(w%A)
    if (allocated(w%v)) deallocate(w%v)
    if (allocated(w%xxt)) deallocate(w%xxt)
    if (allocated(w%extra_scale)) deallocate(w%extra_scale)

    call remove_workspace_evaluate_model(w%evaluate_model_ws,&
         options)

    if (options%model == 4) then 

       call remove_workspace_solve_newton_tensor(& 
            w%solve_newton_tensor_ws, &
            options)
    else

       select case (options%nlls_method)

       case (1) ! use the dogleg method
          call remove_workspace_dogleg(w%dogleg_ws, & 
               options)

       case(2) ! use the AINT method
          call remove_workspace_AINT_tr(w%AINT_tr_ws, & 
               options)

       case(3) ! More-Sorensen 
          call remove_workspace_more_sorensen(&
               w%more_sorensen_ws,options)

       case (4) ! dtrs (Galahad)
          call remove_workspace_solve_galahad(& 
               w%solve_galahad_ws, options)

       end select

    end if

    if (options%scale > 0) then
       call remove_workspace_apply_scaling(w%apply_scaling_ws,options)
    end if

    w%allocated = .false.


  end subroutine remove_workspace_calculate_step

  subroutine setup_workspace_dogleg(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( dogleg_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%d_sd(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%d_gn(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000           
    allocate(w%ghat(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%Jg(m),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    ! setup space for 
    !   solve_LLS
    call setup_workspace_solve_LLS(n,m,w%solve_LLS_ws,options,inform)
    if (inform%status > 0 ) goto 9010
    ! setup space for 
    !   evaluate_model
    call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,options,inform)
    if (inform%status > 0 ) goto 9010

    w%allocated = .true.

    return

    ! Error statements
9000 continue  ! bad allocations in this subroutine
    inform%bad_alloc = 'setup_workspace_dogleg'
    inform%status = ERROR%ALLOCATION
    return

9010 continue  ! bad allocations from dependent subroutine
    return


  end subroutine setup_workspace_dogleg

  subroutine remove_workspace_dogleg(w,options)
    type( dogleg_work ), intent(out) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%d_sd )) deallocate(w%d_sd ) 
    if(allocated( w%d_gn )) deallocate(w%d_gn )
    if(allocated( w%ghat )) deallocate(w%ghat )
    if(allocated( w%Jg )) deallocate(w%Jg )

    ! deallocate space for 
    !   solve_LLS
    call remove_workspace_solve_LLS(w%solve_LLS_ws,options)
    ! deallocate space for 
    !   evaluate_model
    call remove_workspace_evaluate_model(w%evaluate_model_ws,options)

    w%allocated = .false.

    return

  end subroutine remove_workspace_dogleg

  subroutine setup_workspace_solve_LLS(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( solve_LLS_work ) :: w 
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform
    integer :: lwork

    allocate( w%temp(max(m,n)), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    lwork = max(1, min(m,n) + max(min(m,n), 1)*4) 
    allocate( w%work(lwork), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate( w%Jlls(n*m), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    w%allocated = .true.

    return

9000 continue  ! local allocation error
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "solve_LLS"
    return

  end subroutine setup_workspace_solve_LLS

  subroutine remove_workspace_solve_LLS(w,options)
    type( solve_LLS_work ) :: w 
    type( nlls_options ), intent(in) :: options

    if(allocated( w%temp )) deallocate( w%temp)
    if(allocated( w%work )) deallocate( w%work ) 
    if(allocated( w%Jlls )) deallocate( w%Jlls ) 

    w%allocated = .false.

    return

  end subroutine remove_workspace_solve_LLS

  subroutine setup_workspace_evaluate_model(n,m,w,options,inform)
    integer, intent(in) :: n, m        
    type( evaluate_model_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate( w%Jd(m), stat = inform%alloc_status )
    if (inform%alloc_status > 0) goto 9000
    allocate( w%dH(n**2), stat = inform%alloc_status ) 
    if (inform%alloc_status > 0) goto 9000
    allocate( w%dHd(m), stat = inform%alloc_status ) 
    if (inform%alloc_status > 0) goto 9000
    allocate( w%Hd(n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    w%allocated = .true.

    return

9000 continue
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = 'evaluate_model'
    return
  end subroutine setup_workspace_evaluate_model

  subroutine remove_workspace_evaluate_model(w,options)
    type( evaluate_model_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%Jd )) deallocate( w%Jd ) 
    if(allocated( w%dHd )) deallocate( w%dHd ) 
    if(allocated( w%Hd )) deallocate( w%Hd ) 

    w%allocated = .false.

    return

  end subroutine remove_workspace_evaluate_model

  subroutine setup_workspace_AINT_tr(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( AINT_tr_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%B(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%p0(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%p1(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%M0(2*n,2*n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%M1(2*n,2*n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%M0_small(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%M1_small(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%y(2*n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%gtg(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%q(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%LtL(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%y_hardcase(n,2), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    ! setup space for max_eig
    call setup_workspace_max_eig(n,m,w%max_eig_ws,options,inform)
    if (inform%status > 0) goto 9010
    call setup_workspace_evaluate_model(n,m,w%evaluate_model_ws,options,inform)
    if (inform%status > 0) goto 9010
    ! setup space for the solve routine
    if ((options%model .ne. 1)) then
       call setup_workspace_solve_general(n,m,w%solve_general_ws,options,inform)
       if (inform%status > 0 ) goto 9010
    end if

    w%allocated = .true.

    return

9000 continue ! local allocation error
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "AINT_tr"
    !call allocation_error(options,'AINT_tr')
    return

9010 continue ! allocation error from called subroutine
    return

  end subroutine setup_workspace_AINT_tr

  subroutine remove_workspace_AINT_tr(w,options)
    type( AINT_tr_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%B )) deallocate(w%B)
    if(allocated( w%p0 )) deallocate(w%p0)
    if(allocated( w%p1 )) deallocate(w%p1)
    if(allocated( w%M0 )) deallocate(w%M0)
    if(allocated( w%M1 )) deallocate(w%M1)
    if(allocated( w%M0_small )) deallocate(w%M0_small)
    if(allocated( w%M1_small )) deallocate(w%M1_small)
    if(allocated( w%y )) deallocate(w%y)
    if(allocated( w%gtg )) deallocate(w%gtg)
    if(allocated( w%q )) deallocate(w%q)
    if(allocated( w%LtL )) deallocate(w%LtL)
    if(allocated( w%y_hardcase )) deallocate(w%y_hardcase)
    ! setup space for max_eig
    call remove_workspace_max_eig(w%max_eig_ws,options)
    call remove_workspace_evaluate_model(w%evaluate_model_ws,options)
    ! setup space for the solve routine
    if (options%model .ne. 1) then
       call remove_workspace_solve_general(w%solve_general_ws,options)
    end if

    w%allocated = .false.

    return

  end subroutine remove_workspace_AINT_tr

  subroutine setup_workspace_min_eig_symm(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( min_eig_symm_work) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    real(wp), allocatable :: workquery(:)
    integer :: lwork, eigsout

    allocate(w%A(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    allocate(workquery(1),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    if (options%subproblem_eig_fact) then 
       allocate(w%ew(n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000

       call dsyev('V', & ! both ew's and ev's 
            'U', & ! upper triangle of A
            n, w%A, n, & ! data about A
            w%ew, workquery, -1, & 
            inform%external_return)
       if (inform%external_return .ne. 0) goto 9010
    else
       allocate( w%iwork(5*n), stat = inform%alloc_status )
       if (inform%alloc_status > 0) goto 9000
       allocate( w%ifail(n), stat = inform%alloc_status ) 
       if (inform%alloc_status > 0) goto 9000

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
            inform%external_return)
       if (inform%external_return .ne. 0) goto 9020
    end if
    lwork = int(workquery(1))
    if (allocated(workquery)) deallocate(workquery)
    allocate( w%work(lwork), stat = inform%alloc_status )
    if (inform%alloc_status > 0) goto 9000

    w%allocated = .true.

    return

9000 continue      
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "min_eig_symm"
    return

9010 continue
    inform%status = ERROR%FROM_EXTERNAL
    inform%external_name = "lapack_dsyev"

9020 continue
    inform%status = ERROR%FROM_EXTERNAL
    inform%external_name = "lapack_dsyevx"
    return

  end subroutine setup_workspace_min_eig_symm

  subroutine remove_workspace_min_eig_symm(w,options)
    type( min_eig_symm_work) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%A )) deallocate(w%A)        

    if (options%subproblem_eig_fact) then 
       if(allocated( w%ew )) deallocate(w%ew)
    else
       if(allocated( w%iwork )) deallocate( w%iwork )
       if(allocated( w%ifail )) deallocate( w%ifail ) 
    end if
    if(allocated( w%work )) deallocate( w%work ) 

    w%allocated = .false.

    return

  end subroutine remove_workspace_min_eig_symm

  subroutine setup_workspace_max_eig(n,m,w,options,inform)
    integer, intent(in) :: n, m 
    type( max_eig_work) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform), intent(out) :: inform
    real(wp), allocatable :: workquery(:)
    integer :: lwork

    allocate( w%alphaR(2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    w%alphaR = zero
    allocate( w%alphaI(2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    w%alphaI = zero
    allocate( w%beta(2*n),   stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    w%beta = zero
    allocate( w%vr(2*n,2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate( w%ew_array(2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    w%ew_array = zero
    allocate(workquery(1),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    ! make a workspace query to dggev
    call dggev('N', & ! No left eigenvectors
         'V', &! Yes right eigenvectors
         2*n, 1.0, 2*n, 1.0, 2*n, &
         1.0, 0.1, 0.1, & ! eigenvalue data
         0.1, 2*n, & ! not referenced
         0.1, 2*n, & ! right eigenvectors
         workquery, -1, inform%external_return)
    if (inform%external_return > 0) goto 9020
    lwork = int(workquery(1))
    if (allocated(workquery)) deallocate(workquery)
    allocate( w%work(lwork), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate( w%nullindex(2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate( w%vecisreal(2*n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    w%allocated = .true.

    return

9000 continue
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "max_eig"
    return

9020 continue
    inform%status = ERROR%FROM_EXTERNAL
    inform%external_name = "lapack_dggev"
    return

  end subroutine setup_workspace_max_eig

  subroutine remove_workspace_max_eig(w,options)
    type( max_eig_work) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%alphaR )) deallocate( w%alphaR)
    if(allocated( w%alphaI )) deallocate( w%alphaI )
    if(allocated( w%beta )) deallocate( w%beta ) 
    if(allocated( w%vr )) deallocate( w%vr ) 
    if(allocated( w%ew_array )) deallocate( w%ew_array ) 
    if(allocated( w%work )) deallocate( w%work ) 
    if(allocated( w%nullindex )) deallocate( w%nullindex ) 
    if(allocated( w%vecisreal )) deallocate( w%vecisreal )

    w%allocated = .false.

    return

  end subroutine remove_workspace_max_eig

  subroutine setup_workspace_solve_general(n, m, w, options, inform)
    integer, intent(in) :: n, m 
    type( solve_general_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform), intent(out) :: inform

    allocate( w%A(n,n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate( w%ipiv(n), stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    w%allocated = .true.

    return

9000 continue ! allocation error
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "solve_general"
    return

  end subroutine setup_workspace_solve_general

  subroutine remove_workspace_solve_general(w, options)
    type( solve_general_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%A )) deallocate( w%A ) 
    if(allocated( w%ipiv )) deallocate( w%ipiv ) 

    w%allocated = .false.

    return

  end subroutine remove_workspace_solve_general

  subroutine setup_workspace_solve_galahad(n,m,w,options,inform)
    integer, intent(in) :: n,m
    type( solve_galahad_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%ev(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%v_trans(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%ew(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%d_trans(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,options,inform)
    if (inform%status > 0) goto 9010

    w%allocated = .true.

    return

9000 continue ! allocation error here
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "solve_galahad"
    return

9010 continue  ! allocation error from called subroutine
    return

  end subroutine setup_workspace_solve_galahad

  subroutine remove_workspace_solve_galahad(w,options)
    type( solve_galahad_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%ev )) deallocate(w%ev)
    if(allocated( w%v_trans )) deallocate(w%v_trans)
    if(allocated( w%ew )) deallocate(w%ew)
    if(allocated( w%d_trans )) deallocate(w%d_trans)

    call remove_workspace_all_eig_symm(w%all_eig_symm_ws,options)

    w%allocated = .false.

    return

  end subroutine remove_workspace_solve_galahad

  subroutine setup_workspace_all_eig_symm(n,m,w,options,inform)
    integer, intent(in) :: n,m
    type( all_eig_symm_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    real(wp), allocatable :: workquery(:)
    real(wp) :: A,ew
    integer :: lwork

    A = 1.0_wp
    ew = 1.0_wp

    allocate(workquery(1), stat = inform%alloc_status)
    if (inform%alloc_status .ne. 0 ) goto 8000
    call dsyev('V', & ! both ew's and ev's 
         'U', & ! upper triangle of A
         n, A, n, & ! data about A
         ew, workquery, -1, & 
         inform%external_return)  
    if (inform%external_return .ne. 0) goto 9000

    lwork = int(workquery(1))
    if (allocated(workquery)) deallocate(workquery)
    allocate( w%work(lwork), stat = inform%alloc_status )
    if (inform%alloc_status > 0) goto 8000

    w%allocated = .true.

    return

8000 continue  ! allocation error
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "all_eig_sym"
    return

9000 continue ! error from lapack
    inform%status = ERROR%FROM_EXTERNAL
    inform%external_name = "lapack_dsyev"
    return

  end subroutine setup_workspace_all_eig_symm

  subroutine remove_workspace_all_eig_symm(w,options)
    type( all_eig_symm_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%work )) deallocate( w%work ) 

    w%allocated = .false.

  end subroutine remove_workspace_all_eig_symm

  subroutine setup_workspace_more_sorensen(n,m,w,options,inform)
    integer, intent(in) :: n,m
    type( more_sorensen_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%LtL(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%q(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%y1(n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000
    allocate(w%AplusSigma(n,n),stat = inform%alloc_status)
    if (inform%alloc_status > 0) goto 9000

    call setup_workspace_min_eig_symm(n,m,w%min_eig_symm_ws,options,inform)
    if (inform%status > 0) goto 9010

    w%allocated = .true.

    return

9000 continue ! allocation error here
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "more_sorenesen"
    return

9010 continue ! error from called subroutine
    return

  end subroutine setup_workspace_more_sorensen

  subroutine remove_workspace_more_sorensen(w,options)
    type( more_sorensen_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%LtL )) deallocate(w%LtL)
    if(allocated( w%q )) deallocate(w%q)
    if(allocated( w%y1 )) deallocate(w%y1)
    if(allocated( w%AplusSigma )) deallocate(w%AplusSigma)

    call remove_workspace_min_eig_symm(w%min_eig_symm_ws,options)

    w%allocated = .false.

    return

  end subroutine remove_workspace_more_sorensen


  subroutine setup_workspace_apply_scaling(n,m,w,options,inform)

    integer, intent(in) :: n,m
    type( apply_scaling_work ) :: w
    type( nlls_options ), intent(in) :: options
    type( nlls_inform ), intent(out) :: inform

    allocate(w%diag(n), stat = inform%alloc_status )
    if (inform%alloc_status .ne. 0) goto 1000
    w%diag = one
    allocate(w%ev(n,n), stat = inform%alloc_status )
    if (inform%alloc_status .ne. 0) goto 1000

    if (options%scale == 4) then
       allocate(w%tempvec(n))
       call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,options,inform)
       if (inform%status .ne. 0) goto 1010
    end if

    w%allocated = .true.

    return

1000 continue ! allocation error here
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = "apply_scaling"
    return

1010 continue ! error from lower down subroutine
    return

  end subroutine setup_workspace_apply_scaling

  subroutine remove_workspace_apply_scaling(w,options)
    type( apply_scaling_work ) :: w
    type( nlls_options ), intent(in) :: options

    if(allocated( w%diag )) deallocate( w%diag )
    if(allocated( w%ev )) deallocate( w%ev ) 

    w%allocated = .false.

    return 

  end subroutine remove_workspace_apply_scaling

end module ral_nlls_workspaces
