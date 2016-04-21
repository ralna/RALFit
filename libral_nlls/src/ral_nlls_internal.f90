! ral_nlls_internal :: internal subroutines for ral_nlls

module ral_nlls_internal

  use RAL_NLLS_DTRS_double
  use RAL_NLLS_DRQS_double
  
  implicit none

  private

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


  TYPE NLLS_ERROR
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

  END TYPE NLLS_ERROR

  type(NLLS_ERROR) :: ERROR

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
!     norm of the gradient of the objective function is smaller than 
!       MAX( %stop_g_absolute, %stop_g_relative * norm of the initial gradient
!     or if the step is less than %stop_s

     REAL ( KIND = wp ) :: stop_g_absolute = tenm5
     REAL ( KIND = wp ) :: stop_g_relative = tenm8
     
!   should we scale the initial trust region radius?
     
     integer :: relative_tr_radius = 0

!   if relative_tr_radius == 1, then pick a scaling parameter
!   Madsen, Nielsen and Tingleff say pick this to be 1e-6, say, if x_0 is good,
!   otherwise 1e-3 or even 1 would be good starts...
     
     real (kind = wp) :: initial_radius_scale = 1.0

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

     logical :: calculate_svd_J = .true.

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! O U T P U T   C O N T R O L S !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Shall we output progess vectors at termination of the routine?
     logical :: output_progress_vectors = .false.

     

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
  
  type, public :: params_base_type
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
       integer :: nullevs_cols
       real(wp), allocatable :: nullevs(:,:)
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

    type, private :: apply_scaling_work ! workspace for subrouine apply_scaling
       real(wp), allocatable :: diag(:)
       real(wp), allocatable :: ev(:,:)
       real(wp), allocatable :: tempvec(:)
       type( all_eig_symm_work ) :: all_eig_symm_ws
    end type apply_scaling_work
        
    type, private :: solve_galahad_work ! workspace for subroutine dtrs_work
       real(wp), allocatable :: A(:,:), ev(:,:), ew(:), v(:), v_trans(:), d_trans(:)
       type( all_eig_symm_work ) :: all_eig_symm_ws
       type( apply_scaling_work ) :: apply_scaling_ws
    end type solve_galahad_work
        
    type, private :: more_sorensen_work ! workspace for subroutine more_sorensen
 !      type( solve_spd_work ) :: solve_spd_ws
       real(wp), allocatable :: A(:,:), LtL(:,:), AplusSigma(:,:)
       real(wp), allocatable :: v(:), q(:), y1(:)
!       type( solve_general_work ) :: solve_general_ws
       type( min_eig_symm_work ) :: min_eig_symm_ws
       type( apply_scaling_work ) :: apply_scaling_ws
    end type more_sorensen_work

    type, private :: AINT_tr_work ! workspace for subroutine AINT_tr
       type( max_eig_work ) :: max_eig_ws
       type( evaluate_model_work ) :: evaluate_model_ws
       type( solve_general_work ) :: solve_general_ws
!       type( solve_spd_work ) :: solve_spd_ws
       REAL(wp), allocatable :: A(:,:), LtL(:,:), v(:), B(:,:), p0(:), p1(:)
       REAL(wp), allocatable :: M0(:,:), M1(:,:), y(:), gtg(:,:), q(:)
       REAL(wp), allocatable :: M0_small(:,:), M1_small(:,:)
       REAL(wp), allocatable :: y_hardcase(:,:)
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
       type( solve_galahad_work ) :: solve_galahad_ws
    end type calculate_step_work

    type, private :: get_svd_J_work ! workspace for subroutine get_svd_J
       real(wp), allocatable :: Jcopy(:), S(:), work(:)

    end type get_svd_J_work

    type, public :: NLLS_workspace ! all workspaces called from the top level
       integer :: first_call = 1
       integer :: iter = 0 
       real(wp) :: normF0, normJF0, normF, normJF
       real(wp) :: normJFold, normJF_Newton
       real(wp) :: Delta
       real(wp) :: normd
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
       type ( evaluate_model_work ) :: evaluate_model_ws
       type ( get_svd_J_work ) :: get_svd_J_ws
       real(wp) :: tr_nu = 2.0
       integer :: tr_p = 3
    end type NLLS_workspace

    public :: setup_workspaces, solve_galahad, findbeta, mult_j
    public :: mult_jt, solve_spd, solve_general, matmult_inner
    public :: matmult_outer, outer_product, min_eig_symm, max_eig
    public :: remove_workspaces, get_svd_j, calculate_step, evaluate_model
    public :: update_trust_region_radius, apply_second_order_info, rank_one_update
    public :: test_convergence, calculate_rho
    public :: solve_LLS, shift_matrix
    public :: dogleg, more_sorensen, apply_scaling
    public :: ERROR
    
contains

  SUBROUTINE calculate_step(J,f,hf,g,n,m,Delta,d,normd,options,inform,w)

! -------------------------------------------------------
! calculate_step, find the next step in the optimization
! -------------------------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     TYPE( calculate_step_work ) :: w

     select case (options%nlls_method)
        
     case (1) ! Powell's dogleg
        if (options%print_level >= 2) write(options%out,3000) 'dogleg'
        call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%dogleg_ws)
     case (2) ! The AINT method
        if (options%print_level >= 2) write(options%out,3000) 'AINT_TR'
        call AINT_TR(J,f,hf,n,m,Delta,d,normd,options,inform,w%AINT_tr_ws)
     case (3) ! More-Sorensen
        if (options%print_level >= 2) write(options%out,3000) 'More-Sorensen'
        call more_sorensen(J,f,hf,n,m,Delta,d,normd,options,inform,w%more_sorensen_ws)
     case (4) ! Galahad
        if (options%print_level >= 2) write(options%out,3000) 'DTRS'
        call solve_galahad(J,f,hf,n,m,Delta,d,normd,options,inform,w%solve_galahad_ws)
     case default
        inform%status = ERROR%UNSUPPORTED_METHOD
     end select

     if (options%print_level >= 2) write(options%out,3010)
     
3000 FORMAT('*** Solving the trust region subproblem using ',A,' ***')
3010 FORMAT('*** Trust region subproblem solution found ***')
     

   END SUBROUTINE calculate_step
   
   subroutine apply_scaling(J,n,m,A,v,w,options,inform)
     !-------------------------------
     ! apply_scaling
     ! input :: Jacobian matrix, J
     ! ouput :: scaled Hessisan, H, and J^Tf, v.
     !
     ! Calculates a diagonal scaling W, stored in w%diag
     ! updates v(i) -> (1/W_i) * v(i)
     !         A(i,j) -> (1 / (W_i * W_j)) * A(i,j)
     !-------------------------------
     real(wp), intent(in) :: J(*) 
     integer, intent(in) :: n,m
     real(wp), intent(inout) :: A(:,:)
     real(wp), intent(inout) :: v(:)
     type( apply_scaling_work ), intent(inout) :: w
     type( nlls_options ), intent(in) :: options
     type( nlls_inform ), intent(inout) :: inform

     integer :: ii, jj
     real(wp) :: Jij, temp


     select case (options%scale)
     case (1,2)
        do ii = 1,n
           temp = zero
           what_scale: if (options%scale == 1) then
              ! use the scaling present in gsl:
              ! scale by W, W_ii = ||J(i,:)||_2^2
              do jj = 1,m
                 call get_element_of_matrix(J,m,jj,ii,Jij)
                 temp = temp + Jij**2
              end do
           elseif ( options%scale == 2) then 
              ! scale using the (approximate) hessian
              do jj = 1,n
                 temp = temp + A(ii,jj)**2
              end do
           end if what_scale
           trim_scale: if (temp < options%scale_min) then 
              if (options%scale_trim_min) then 
                 temp = options%scale_min
              else
                 temp = one
              end if
           elseif (temp > options%scale_max) then
              if (options%scale_trim_max) then 
                 temp = options%scale_max
              else
                 temp = one 
              end if
           end if trim_scale
           temp = sqrt(temp)
           if (options%scale_require_increase) then
              w%diag(ii) = max(temp,w%diag(ii))
           else
              w%diag(ii) = temp
           end if
        end do
!!$     case (3)
!!$        ! assuming problems are small...do an eigen-decomposition of H
!!$        write(*,*) '***** Warning ********'
!!$        write(*,*) '*    not robust      *'
!!$        write(*,*) '**********************'
!!$        call all_eig_symm(A,n,w%tempvec,w%ev,w%all_eig_symm_ws,inform)
!!$        if (inform%status .ne. 0) goto 1000
!!$        do ii = 1,n
!!$           ! plain version...
!!$           w%diag(n + 1 - ii) = w%tempvec(ii)
!!$        end do
!!$        ! todo : require_increase, test for trimming
     case default
        inform%status = ERROR%BAD_SCALING
        return
     end select
          
     ! now we have the w%diagonal scaling matrix, actually scale the 
     ! Hessian approximation and J^Tf
     do ii = 1,n
        v(ii) = v(ii) / w%diag(ii)
        A(ii,:) = A(ii,:) / w%diag(ii)
        A(:,ii) = A(:,ii) / w%diag(ii)
     end do
   
     return
     
1000 continue
     ! error in external package
     return

   end subroutine apply_scaling

   SUBROUTINE dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w)
! -----------------------------------------
! dogleg, implement Powell's dogleg method
! -----------------------------------------

     REAL(wp), intent(in) :: J(:), hf(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     TYPE( dogleg_work ) :: w
     
     real(wp) :: alpha, beta

     !     Jg = J * g
     call mult_J(J,n,m,g,w%Jg)

     alpha = norm2(g)**2 / norm2( w%Jg )**2
       
     w%d_sd = alpha * g;

     ! Solve the linear problem...
     select case (options%model)
     case (1)
        ! linear model...
        call solve_LLS(J,f,n,m,w%d_gn,inform,w%solve_LLS_ws)
        if ( inform%status .ne. 0 ) goto 1000
     case default
        inform%status = ERROR%DOGLEG_MODEL
        return
     end select
     
     if (norm2(w%d_gn) <= Delta) then
        d = w%d_gn
        if (options%print_level >=2) write(options%out,2000)
     else if (norm2( alpha * w%d_sd ) >= Delta) then
        d = (Delta / norm2(w%d_sd) ) * w%d_sd
        if (options%print_level >=2) write(options%out,2010)
     else
        w%d_sd = alpha * w%d_sd
        w%ghat = w%d_gn - w%d_sd
        call findbeta(w%d_sd,w%ghat,Delta,beta,inform)
        if ( inform%status .ne. 0 ) goto 1000
        d = w%d_sd + beta * w%ghat
        if (options%print_level >=2) write(options%out,2020)
     end if

     normd = norm2(d)
     
     return
     
1000 continue 
     ! bad error return from solve_LLS
     return

! Printing commands
2000 FORMAT('Gauss Newton step taken')
2010 FORMAT('Steepest descent step taken')
2020 FORMAT('Dogleg step taken')

   END SUBROUTINE dogleg
     
   SUBROUTINE AINT_tr(J,f,hf,n,m,Delta,d,normd,options,inform,w)
     ! -----------------------------------------
     ! AINT_tr
     ! Solve the trust-region subproblem using 
     ! the method of ADACHI, IWATA, NAKATSUKASA and TAKEDA
     ! -----------------------------------------

     REAL(wp), intent(in) :: J(:), f(:), hf(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( AINT_tr_work ) :: w
        
     integer :: keep_p0, i, size_hard(2)
     real(wp) :: obj_p0, obj_p1
     REAL(wp) :: norm_p0, tau, lam, eta
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
     
     select case (options%model)
     case (1)
        call solve_spd(w%A,-w%v,w%LtL,w%p0,n,inform)
        if (inform%status .ne. 0) goto 1000
     case default
        call solve_general(w%A,-w%v,w%p0,n,inform,w%solve_general_ws)
        if (inform%status .ne. 0) goto 1000
     end select
          
     call matrix_norm(w%p0,w%B,norm_p0)
     
     if (norm_p0 < Delta) then
        keep_p0 = 1;
        ! get obj_p0 : the value of the model at p0
        if (options%print_level >=3) write(options%out,2000) 'p0'     
        call evaluate_model(f,J,hf,w%p0,obj_p0,m,n,options,w%evaluate_model_ws)
     end if

     w%M0(1:n,1:n) = -w%B
     w%M0(n+1:2*n,1:n) = w%A
     w%M0(1:n,n+1:2*n) = w%A
     call outer_product(w%v,n,w%gtg) ! gtg = Jtg * Jtg^T
     w%M0(n+1:2*n,n+1:2*n) = (-1.0 / Delta**2) * w%gtg

     w%M1 = 0.0
     w%M1(n+1:2*n,1:n) = -w%B
     w%M1(1:n,n+1:2*n) = -w%B
     
     call max_eig(w%M0,w%M1,2*n,lam, w%y, w%y_hardcase, options, inform, w%max_eig_ws)
     if ( inform%status > 0 ) goto 1000

     if (norm2(w%y(1:n)) < tau) then
        ! Hard case
        if ( options%print_level >= 3) write(options%out, 2010)
        ! overwrite H onto M0, and the outer prod onto M1...
        size_hard = shape(w%y_hardcase)
        call matmult_outer( matmul(w%B,w%y_hardcase), size_hard(2), n, w%M1_small)
        w%M0_small = w%A(:,:) + lam*w%B(:,:) + w%M1_small
        ! solve Hq + g = 0 for q
        select case (options%model) 
        case (1)
           call solve_spd(w%M0_small,-w%v,w%LtL,w%q,n,inform)
           if (inform%status .ne. 0) goto 1000
        case default
          call solve_general(w%M0_small,-w%v,w%q,n,inform,w%solve_general_ws)
          if (inform%status .ne. 0) goto 1000
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! (I think..) and inside...fix

        
        ! find max eta st ||q + eta v(:,1)||_B = Delta
        call findbeta(w%q,w%y_hardcase(:,1),Delta,eta,inform)
        if ( inform%status .ne. 0 ) goto 1000

        !!!!!      ^^TODO^^    !!!!!
        ! currently assumes B = I !!
        !!!!       fixme!!      !!!!
        
        w%p1(:) = w%q(:) + eta * w%y_hardcase(:,1)
        
     else 
        select case (options%model)
        case (1)
           call solve_spd(w%A + lam*w%B,-w%v,w%LtL,w%p1,n,inform)
           if (inform%status .ne. 0) goto 1000
        case default
           call solve_general(w%A + lam*w%B,-w%v,w%p1,n,inform,w%solve_general_ws)
           if (inform%status .ne. 0) goto 1000
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! and inside...fix
     end if
     
     ! get obj_p1 : the value of the model at p1
     if (options%print_level >=3) write(options%out,2000) 'p1'     
     call evaluate_model(f,J,hf,w%p1,obj_p1,m,n,options,w%evaluate_model_ws)

     ! what gives the smallest objective: p0 or p1?
     if (obj_p0 < obj_p1) then
        d = w%p0
        if (options%print_level >=2) write(options%out,2030) 'p0'
     else 
        d = w%p1
        if (options%print_level >=2) write(options%out,2030) 'p1'
     end if

     normd = norm2(d)

     return
         
1000 continue 
     ! bad error return from external package
     return

! print statements   
2000 FORMAT('Evaluating the model at ',A2,':')
2010 FORMAT('Hard case identified')
2030 FORMAT(A2,' chosen as d')
 

   END SUBROUTINE AINT_tr

   subroutine more_sorensen(J,f,hf,n,m,Delta,d,nd,options,inform,w)
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
     real(wp), intent(out) :: nd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ) :: w

     real(wp) :: nq, epsilon
     real(wp) :: sigma, alpha, local_ms_shift, sigma_shift
     integer :: i, no_restarts

     
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

     ! if scaling needed, do it
     if ( options%scale .ne. 0) then
        call apply_scaling(J,n,m,w%A,w%v,w%apply_scaling_ws,options,inform)
     end if
     
     local_ms_shift = options%more_sorensen_shift

     ! d = -A\v
     call solve_spd(w%A,-w%v,w%LtL,d,n,inform)
     if (inform%status .eq. 0) then
        ! A is symmetric positive definite....
        sigma = zero
        if (options%print_level >=3) write(options%out,6000)
     else
        ! reset the error calls -- handled in the code....
        inform%status = 0
        inform%external_return = 0
        inform%external_name = REPEAT( ' ', 80 )
        call min_eig_symm(w%A,n,sigma,w%y1,options,inform,w%min_eig_symm_ws) 
        if (inform%status .ne. 0) goto 1000
        sigma = -(sigma - local_ms_shift)
        if (options%print_level >= 3) write(options%out,6010) sigma
        ! find a shift that makes (A + sigma I) positive definite
        call get_pd_shift(n,sigma,d,options,inform,w)
        if (inform%status .ne. 0) goto 4000
        if (options%print_level >=3) write(options%out,6020)
     end if
     
     nd = norm2(d)
     
     if (options%print_level >= 2) write(options%out,5000)
     ! now, we're not in the trust region initally, so iterate....
     sigma_shift = zero
     no_restarts = 0
     ! set 'small' in the context of the algorithm
     epsilon = max( options%more_sorensen_tol * Delta, options%more_sorensen_tiny )
     do i = 1, options%more_sorensen_maxits
        if (options%print_level >= 2) write(options%out,5010) i-1, nd, sigma, sigma_shift
                
        if (nd .le. Delta + epsilon) then
           ! we're within the tr radius
           if (options%print_level >= 3) write(options%out,6030)
           if ( abs(sigma) < options%more_sorensen_tiny ) then
              ! we're good....exit
              if (options%print_level >= 3) write(options%out,6040)
              goto 1020
           else if ( abs( nd - Delta ) < epsilon ) then
              ! also good...exit
              if (options%print_level >= 3) write(options%out,6050)
              goto 1020              
           end if
           call findbeta(d,w%y1,Delta,alpha,inform)
           if (inform%status .ne. 0 ) goto 1000  
           d = d + alpha * w%y1
           if (options%print_level >= 3) write(options%out,6060)
           ! also good....exit
           goto 1020
        end if

        w%q = d ! w%q = R'\d
        CALL DTRSM( 'Left', 'Lower', 'No Transpose', 'Non-unit', n, & 
             1, one, w%LtL, n, w%q, n )
        
        nq = norm2(w%q)
        if (options%print_level >= 3) write(options%out,6080) nq
        
        sigma_shift = ( (nd/nq)**2 ) * ( (nd - Delta) / Delta )
        if (abs(sigma_shift) < options%more_sorensen_tiny * abs(sigma) ) then
           if (no_restarts < 1) then 
              ! find a shift that makes (A + sigma I) positive definite
              call get_pd_shift(n,sigma,d,options,inform,w)
              if (inform%status .ne. 0) goto 4000
              no_restarts = no_restarts + 1
           else
              ! we're not going to make progress...jump out 
              inform%status = ERROR%MS_NO_PROGRESS
              goto 4000
           end if
        else 
           sigma = sigma + sigma_shift
        end if

        call shift_matrix(w%A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-w%v,w%LtL,d,n,inform)
        if (inform%status .ne. 0) goto 1000
        
        nd = norm2(d)

     end do
     if (options%print_level >= 2) write(options%out,5010)

     goto 1040
     
1000 continue 
     ! bad error return from external package
     goto 4000
     
1020 continue
     ! inital point was successful
     if (options%print_level==2) write(options%out,5040)
     goto 4000

1040 continue
     ! maxits reached, not converged
     if (options%print_level >=2) write(options%out,5020)
     inform%status = ERROR%MS_MAXITS
     goto 4000

3000 continue
     ! too many shifts
     inform%status = ERROR%MS_TOO_MANY_SHIFTS
     goto 4000
     
4000 continue
     ! exit the routine
     if (options%scale .ne. 0 ) then 
        do i = 1, n
           d(i) = d(i) / w%apply_scaling_ws%diag(i)
        end do
     end if
     return 

! Printing statements
! print_level >= 2 
5000 FORMAT('iter',4x,'nd',12x,'sigma',9x,'sigma_shift')
5010 FORMAT(i4,2x,ES12.4,2x,ES12.4,2x,ES12.4)
5020 FORMAT('More-Sorensen failed to converge within max number of iterations')   
5030 FORMAT('More-Sorensen converged at iteration ',i4)
5040 FORMAT('Leaving More-Sorensen')
  
! print_level >= 3 
6000 FORMAT('A is symmetric positive definite')     
6010 FORMAT('Trying a shift of sigma = ',ES12.4)     
6020 FORMAT('A + sigma I is symmetric positive definite') 
6030 FORMAT('We''re within the trust region radius initially')     
6040 FORMAT('Sigma tiny, so exit')  
6050 FORMAT('||d|| = Delta, so exit')   
6060 FORMAT('Return d + alpha*y_1') 
6070 FORMAT('Converged! ||d|| = Delta at iteration ',i4)      
6080 FORMAT('nq = ',ES12.4)

   end subroutine more_sorensen

   subroutine get_pd_shift(n,sigma,d,options,inform,w)

     !--------------------------------------------------
     ! get_pd_shift
     !
     ! Given an indefinite matrix w%A, find a shift sigma
     ! such that (A + sigma I) is positive definite
     !--------------------------------------------------
     
     integer, intent(in) :: n 
     real(wp), intent(inout) :: sigma
     real(wp), intent(inout) :: d(:)
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ), intent(inout) :: w

     integer :: no_shifts
     logical :: successful_shift
     
     no_shifts = 0
     successful_shift = .false.
     do while( .not. successful_shift )
        call shift_matrix(w%A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-w%v,w%LtL,d,n,inform)
        if ( inform%status .ne. 0 ) then
           ! reset the error calls -- handled in the code....
           inform%status = 0
           inform%external_return = 0
           inform%external_name = REPEAT( ' ', 80 )
           no_shifts = no_shifts + 1
           if ( no_shifts == 10 ) goto 3000 ! too many shifts -- exit
           sigma =  sigma + (10**no_shifts) * options%more_sorensen_shift
           if (options%print_level >=3) write(options%out,6010) sigma
        else
           successful_shift = .true.
        end if
     end do

     return

3000 continue
     ! too many shifts
     inform%status = ERROR%MS_TOO_MANY_SHIFTS
     return     

6010 FORMAT('Trying a shift of sigma = ',ES12.4)
     

   end subroutine get_pd_shift
   
   subroutine solve_galahad(J,f,hf,n,m,Delta,d,normd,options,inform,w)

     !---------------------------------------------
     ! solve_galahad
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
     real(wp), intent(out) :: normd ! ||d||_D, where D is the scaling
     type( solve_galahad_work ) :: w
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform

     TYPE ( DTRS_CONTROL_TYPE ) :: dtrs_options
     TYPE ( DTRS_inform_type )  :: dtrs_inform
     TYPE ( DRQS_CONTROL_TYPE ) :: drqs_options
     TYPE ( DRQS_inform_type )  :: drqs_inform

!     real(wp), allocatable :: diag(:)
     integer :: ii

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

     ! if scaling needed, do it
     if ( options%scale .ne. 0) then
        call apply_scaling(J,n,m,w%A,w%v,w%apply_scaling_ws,options,inform)
     end if

     ! Now that we have the unprocessed matrices, we need to get an 
     ! eigendecomposition to make A diagonal
     !
     call all_eig_symm(w%A,n,w%ew,w%ev,w%all_eig_symm_ws,inform)
     if (inform%status .ne. 0) goto 1000

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
      
     do ii = 1,n
        if (abs(w%v_trans(ii)) < epsmch) then
           w%v_trans(ii) = zero
        end if
        if (abs(w%ew(ii)) < epsmch) then
           w%ew(ii) = zero
        end if
     end do
 
    select case (options%type_of_method)
     case (1)
        call dtrs_initialize( dtrs_options, dtrs_inform ) 
        dtrs_options%error = options%error
        dtrs_options%out = options%out
        dtrs_options%print_level = options%print_level - 1
        call dtrs_solve(n, Delta, zero, w%v_trans, w%ew, w%d_trans, dtrs_options, dtrs_inform )
        if ( dtrs_inform%status .ne. 0) then
           inform%external_return = dtrs_inform%status
           inform%external_name = 'galahad_dtrs'
           inform%status = ERROR%FROM_EXTERNAL
           goto 1000
        end if
     case(2)
        call drqs_initialize( drqs_options, drqs_inform ) 
        drqs_options%error = options%error
        drqs_options%out = options%out
        drqs_options%print_level = options%print_level - 1
        call drqs_solve & 
             (n,2.0_wp,1/Delta, zero, w%v_trans, w%ew, w%d_trans, drqs_options, drqs_inform)
        if ( drqs_inform%status .ne. 0) then
           inform%external_return = drqs_inform%status
           inform%external_name = 'galahad_drqs'
           inform%status = ERROR%FROM_EXTERNAL
           goto 1000
        end if
     end select
  ! and return the un-transformed vector
     call mult_J(w%ev,n,n,w%d_trans,d)

     normd = norm2(d) ! ||d||_D
     
     if (options%scale .ne. 0 ) then 
        do ii = 1, n
           d(ii) = d(ii) / w%apply_scaling_ws%diag(ii)
        end do
     end if

     return

1000 continue 
     ! bad error return from external package
     return
     
   end subroutine solve_galahad


   SUBROUTINE solve_LLS(J,f,n,m,d_gn,inform,w)
       
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------

       REAL(wp), DIMENSION(:), INTENT(IN) :: J
       REAL(wp), DIMENSION(:), INTENT(IN) :: f
       INTEGER, INTENT(IN) :: n, m
       REAL(wp), DIMENSION(:), INTENT(OUT) :: d_gn
       type(NLLS_inform), INTENT(INOUT) :: inform

       character(1) :: trans = 'N'
       integer :: nrhs = 1, lwork, lda, ldb
       type( solve_LLS_work ) :: w
       
       
       lda = m
       ldb = max(m,n)
       w%temp(1:m) = f(1:m)
       lwork = size(w%work)
       
       w%Jlls(:) = J(:)
       
       call dgels(trans, m, n, nrhs, w%Jlls, lda, w%temp, ldb, w%work, lwork, &
            inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dgels'
          return
       end if

       d_gn = -w%temp(1:n)
              
     END SUBROUTINE solve_LLS
     
     SUBROUTINE findbeta(a, b, Delta, beta, inform)

!  -----------------------------------------------------------------
!  findbeta, a subroutine to find the optimal beta such that 
!   || d || = Delta, where d = a + beta * b
!   
!   uses the approach from equation (3.20b), 
!    "Methods for non-linear least squares problems" (2nd edition, 2004)
!    by Madsen, Nielsen and Tingleff      
!  -----------------------------------------------------------------

     real(wp), dimension(:), intent(in) :: a, b 
     real(wp), intent(in) ::  Delta
     real(wp), intent(out) :: beta
     type( nlls_inform ), intent(out) :: inform
     
     real(wp) :: c, normb2, norma2, discrim, denom
     
     c = dot_product(a,b)

     norma2 = norm2(a)**2
     normb2 = norm2(b)**2

     discrim = c**2 + (normb2)*(Delta**2 - norma2);
     if ( discrim < zero ) then
        inform%status = ERROR%FIND_BETA
        inform%external_name = 'findbeta'
        return
     end if

     if (c .le. 0) then
        beta = (-c + sqrt(discrim) ) / normb2
     else
        beta = (Delta**2 - norma2) / ( c + sqrt(discrim) )
     end if
        

     END SUBROUTINE findbeta
     
     subroutine evaluate_model(f,J,hf,d,md,m,n,options,w)
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
       TYPE( nlls_options ), INTENT( IN ) :: options
       type( evaluate_model_work ) :: w
       
       !Jd = J*d
       call mult_J(J,n,m,d,w%Jd)
       
       ! First, get the base 
       ! 0.5 (f^T f + f^T J d + d^T' J ^T J d )
       md = 0.5 * norm2(f + w%Jd)**2
       select case (options%model)
       case (1) ! first-order (no Hessian)
          ! nothing to do here...
          continue
       case default
          ! these have a dynamic H -- recalculate
          ! H = J^T J + HF, HF is (an approx?) to the Hessian
          call mult_J(hf,n,n,d,w%Hd)
          md = md + 0.5 * dot_product(d,w%Hd)
       end select
       if (options%print_level >= 3) write(options%out,1000) md

1000   FORMAT('Model evauated successfully.  m_k(d) = ',ES12.4)

     end subroutine evaluate_model
     
     subroutine calculate_rho(normf,normfnew,md,rho,options)
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
       TYPE( nlls_options ), INTENT( IN ) :: options

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

       if (options%print_level >= 3) write(options%out,1000) actual_reduction
       if (options%print_level >= 3) write(options%out,1010) predicted_reduction
       if (options%print_level >= 3) write(options%out,1020) rho
       
1000   FORMAT('Actual reduction (in cost function) = ', ES12.4)
1010   FORMAT('Predicted reduction (in model) = ', ES12.4)
1020   FORMAT('rho returned = ', ES12.4)       

     end subroutine calculate_rho

     subroutine apply_second_order_info(n,m,&
          X,w, & !f,hf,
          eval_Hf,&
          !          d, y, y_sharp, & 
          params,options,inform,weights)
       integer, intent(in)  :: n, m 
       real(wp), intent(in) :: X(:)!, f(:)
!       real(wp), intent(inout) :: hf(:)
!       real(wp), intent(in) :: d(:), y(:), y_sharp(:)
       type( NLLS_workspace ), intent(inout) :: w
       procedure( eval_hf_type ) :: eval_Hf
       class( params_base_type ) :: params
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       real(wp), intent(in), optional :: weights(:)
       
       if (options%exact_second_derivatives) then
          if ( present(weights) ) then
             call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
          else
             call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
          end if
          inform%h_eval = inform%h_eval + 1
       else
          ! use the rank-one approximation...
          call rank_one_update(w%hf,w,n)                      
       end if

     end subroutine apply_second_order_info


     subroutine rank_one_update(hf,w,n)

       real(wp), intent(inout) :: hf(:)
       type( NLLS_workspace ), intent(inout) :: w
       integer, intent(in) :: n
      
       real(wp) :: yts, alpha, dSks

       yts = dot_product(w%d,w%y)
       if ( abs(yts) < 10 * epsmch ) then
          ! safeguard: skip this update
          return
       end if

       call mult_J(hf,n,n,w%d,w%Sks) ! hfs = S_k * d

       w%ysharpSks = w%y_sharp - w%Sks

       ! now, let's scale hd (Nocedal and Wright, Section 10.2)
       dSks = abs(dot_product(w%d,w%Sks))
       alpha = abs(dot_product(w%d,w%y_sharp))/ dSks
       alpha = min(one,alpha)
       hf(:)  = alpha * hf(:)

       ! update S_k (again, as in N&W, Section 10.2)

       ! hf = hf + (1/yts) (y# - Sk d)^T y:
       alpha = 1/yts
       call dGER(n,n,alpha,w%ysharpSks,1,w%y,1,hf,n)
       ! hf = hf + (1/yts) y^T (y# - Sk d):
       call dGER(n,n,alpha,w%y,1,w%ysharpSks,1,hf,n)
       ! hf = hf - ((y# - Sk d)^T d)/((yts)**2)) * y y^T
       alpha = -dot_product(w%ysharpSks,w%d)/(yts**2)
       call dGER(n,n,alpha,w%y,1,w%y,1,hf,n)

     end subroutine rank_one_update


     subroutine update_trust_region_radius(rho,options,inform,w)

       real(wp), intent(inout) :: rho ! ratio of actual to predicted reduction
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_workspace ), intent(inout) :: w

       select case(options%tr_update_strategy)
       case(1) ! default, step-function
          if (rho < options%eta_success_but_reduce) then
             ! unsuccessful....reduce Delta
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             if (options%print_level > 2) write(options%out,3010) w%Delta     
          else if (rho < options%eta_very_successful) then 
             ! doing ok...retain status quo
             if (options%print_level > 2) write(options%out,3020) w%Delta 
          else if (rho < options%eta_too_successful ) then
             ! more than very successful -- increase delta
             w%Delta = min(options%maximum_radius, &
                  options%radius_increase * w%normd)
             ! increase based on normd = ||d||_D
             ! if d is on the tr boundary, this is Delta
             ! otherwise, point was within the tr, and there's no point increasing 
             ! the radius
             if (options%print_level > 2) write(options%out,3030) w%Delta
          else if (rho >= options%eta_too_successful) then
             ! too successful....accept step, but don't change w%Delta
             if (options%print_level > 2) write(options%out,3040) w%Delta 
          else
             ! just incase (NaNs and the like...)
             if (options%print_level > 2) write(options%out,3050) w%Delta 
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             rho = -one ! set to be negative, so that the logic works....
          end if
       case(2) ! Continuous method
          ! Based on that proposed by Hans Bruun Nielsen, TR IMM-REP-1999-05
          ! http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf
          if (rho >= options%eta_too_successful) then
             ! too successful....accept step, but don't change w%Delta
             if (options%print_level > 2) write(options%out,3040) w%Delta 
          else if (rho > options%eta_successful) then 
             w%Delta = w%Delta * min(options%radius_increase, &
                  max(options%radius_reduce, & 
                  1 - ( (options%radius_increase - 1) * ((1 - 2*rho)**w%tr_p)) ))
             w%tr_nu = options%radius_reduce
             if (options%print_level > 2) write(options%out,3060) w%Delta 
          else if ( rho <= options%eta_successful ) then 
             w%Delta = w%Delta * w%tr_nu
             w%tr_nu =  w%tr_nu * 0.5_wp
             if (options%print_level > 2) write(options%out,3010) w%Delta
          else
             ! just incase (NaNs and the like...)
             if (options%print_level > 2) write(options%out,3050) w%Delta 
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta 
             rho = -one ! set to be negative, so that the logic works....
          end if
       case default
          inform%status = ERROR%BAD_TR_STRATEGY
          return          
       end select

       return

       ! print statements

3010   FORMAT('Unsuccessful step -- decreasing Delta to', ES12.4)      
3020   FORMAT('Successful step -- Delta staying at', ES12.4)     
3030   FORMAT('Very successful step -- increasing Delta to', ES12.4)
3040   FORMAT('Step too successful -- Delta staying at', ES12.4) 
3050   FORMAT('NaN encountered -- reduced Delta to', ES12.4)   
3060   FORMAT('Changing Delta to ', ES12.4)


     end subroutine update_trust_region_radius

     subroutine test_convergence(normF,normJF,normF0,normJF0,options,inform)

       real(wp), intent(in) :: normF, normJf, normF0, normJF0
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform

       if ( normF <= max(options%stop_g_absolute, &
            options%stop_g_relative * normF0) ) then
          inform%convergence_normf = 1
          return
       end if

       if ( (normJF/normF) <= max(options%stop_g_absolute, &
            options%stop_g_relative * (normJF0/normF0)) ) then
          inform%convergence_normg = 1
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

     subroutine get_element_of_matrix(J,m,ii,jj,Jij)
       real(wp), intent(in) :: J(*)
       integer, intent(in) :: m
       integer, intent(in) :: ii,jj
       real(wp), intent(out) :: Jij

       ! return the (ii,jj)th entry of a matrix 

       ! J held by columns....
       Jij = J(ii + (jj-1)*m)

     end subroutine get_element_of_matrix

     subroutine solve_spd(A,b,LtL,x,n,inform)
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: LtL(:,:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform), intent(inout) :: inform

       ! A wrapper for the lapack subroutine dposv.f
       ! get workspace for the factors....
       LtL(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dposv('L', n, 1, LtL, n, x, n, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dposv'
          return
       end if

     end subroutine solve_spd

     subroutine solve_general(A,b,x,n,inform,w)
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform ), intent(inout) :: inform
       type( solve_general_work ) :: w

       ! A wrapper for the lapack subroutine dposv.f
       ! NOTE: A would be destroyed
       w%A(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dgesv( n, 1, w%A, n, w%ipiv, x, n, inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dgesv'
          return
       end if

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

     subroutine all_eig_symm(A,n,ew,ev,w,inform)
       ! calculate all the eigenvalues of A (symmetric)

       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew(:), ev(:,:)
       type( all_eig_symm_work ) :: w
       type( nlls_inform ), intent(inout) :: inform

       integer :: lwork

       ! copy the matrix A into the eigenvector array
       ev(1:n,1:n) = A(1:n,1:n)

       lwork = size(w%work)
       ! call dsyev --> all eigs of a symmetric matrix

       call dsyev('V', & ! both ew's and ev's 
            'U', & ! upper triangle of A
            n, ev, n, & ! data about A
            ew, w%work, lwork, & 
            inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dsyev'
          return
       end if

     end subroutine all_eig_symm

     subroutine min_eig_symm(A,n,ew,ev,options,inform,w)
       ! calculate the leftmost eigenvalue of A

       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew, ev(:)
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_options ), INTENT( IN ) :: options
       type( min_eig_symm_work ) :: w

       real(wp) :: tol, dlamch
       integer :: lwork, eigsout, minindex(1)

       tol = 2*dlamch('S')!1e-15

       w%A(1:n,1:n) = A(1:n,1:n) ! copy A, as workspace for dsyev(x)
       ! note that dsyevx (but not dsyev) only destroys the lower (or upper) part of A
       ! so we could possibly reduce memory use here...leaving for 
       ! ease of understanding for now.

       lwork = size(w%work)
       if ( options%subproblem_eig_fact ) then
          ! call dsyev --> all eigs of a symmetric matrix
          call dsyev('V', & ! both ew's and ev's 
               'U', & ! upper triangle of A
               n, w%A, n, & ! data about A
               w%ew, w%work, lwork, & 
               inform%external_return)
          if (inform%external_return .ne. 0) then
             inform%status = ERROR%FROM_EXTERNAL
             inform%external_name = 'lapack_dsyev'
          end if
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
               inform%external_return)
          if (inform%external_return .ne. 0) then
             inform%status = ERROR%FROM_EXTERNAL
             inform%external_name = 'lapack_dsyevx'
          end if
       end if

       return

     end subroutine min_eig_symm

     subroutine max_eig(A,B,n,ew,ev,nullevs,options,inform,w)

       real(wp), intent(inout) :: A(:,:), B(:,:)
       integer, intent(in) :: n 
       real(wp), intent(out) :: ew, ev(:)
       real(wp), intent(inout), allocatable :: nullevs(:,:)
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( max_eig_work ) :: w

       integer :: lwork, maxindex(1), no_null, halfn
       real(wp):: tau
       integer :: i 

       ! Find the max eigenvalue/vector of the generalized eigenproblem
       !     A * y = lam * B * y
       ! further, if ||y(1:n/2)|| \approx 0, find and return the 
       ! eigenvectors y(n/2+1:n) associated with this

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
            w%work, lwork, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dggev'
          return
       end if

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
!          allocate(nullevs(halfn,no_null))
          if (no_null > size(nullevs,2)) then
             ! increase the size of the allocated array only if we need to
             if(allocated(nullevs)) deallocate( nullevs )
             allocate( nullevs(halfn,no_null) , stat = inform%alloc_status)
             if (inform%alloc_status > 0) goto 2000
          end if 
          nullevs(1:halfn,1:no_null) = w%vr(halfn+1 : n,w%nullindex(1:no_null))
       end if

       ew = w%alphaR(maxindex(1))/w%beta(maxindex(1))
       ev(:) = w%vr(:,maxindex(1))

       return 

1000   continue 
       inform%status = ERROR%AINT_EIG_IMAG ! Eigs imaginary error
       return

1010   continue
       inform%status = ERROR%AINT_EIG_ODD
       return
       
2000   continue
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "max_eig"
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

     subroutine get_svd_J(n,m,J,s1,sn,options,status,w)
       integer, intent(in) :: n,m 
       real(wp), intent(in) :: J(:)
       real(wp), intent(out) :: s1, sn
       type( nlls_options ) :: options
       integer, intent(out) :: status
       type( get_svd_J_work ) :: w

       !  Given an (m x n)  matrix J held by columns as a vector,
       !  this routine returns the largest and smallest singular values
       !  of J.

       character :: jobu(1), jobvt(1)
       integer :: lwork

       w%Jcopy(:) = J(:)

       jobu  = 'N' ! calculate no left singular vectors
       jobvt = 'N' ! calculate no right singular vectors

       lwork = size(w%work)

       call dgesvd( JOBU, JOBVT, n, m, w%Jcopy, n, w%S, w%S, 1, w%S, 1, & 
            w%work, lwork, status )
       if ( (status .ne. 0) .and. (options%print_level > 3) ) then 
          write(options%error,'(a,i0)') 'Error when calculating svd, dgesvd returned', &
               status
          s1 = -1.0
          sn = -1.0
          ! allow to continue, but warn user and return zero singular values
       else
          s1 = w%S(1)
          sn = w%S(n)
          if (options%print_level > 2) then 
             write(options%out,'(a,es12.4,a,es12.4)') 's1 = ', s1, '    sn = ', sn
             write(options%out,'(a,es12.4)') 'k(J) = ', s1/sn
          end if
       end if

     end subroutine get_svd_J


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!                                                       !!
     !! W O R K S P A C E   S E T U P   S U B R O U T I N E S !!
     !!                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine setup_workspaces(workspace,n,m,options,inform)

       type( NLLS_workspace ), intent(out) :: workspace
       type( nlls_options ), intent(in) :: options
       integer, intent(in) :: n,m
       type( NLLS_inform ), intent(out) :: inform

       workspace%first_call = 0

       workspace%tr_nu = options%radius_increase
       workspace%tr_p = 7

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


       select case (options%nlls_method)

       case (1) ! use the dogleg method
          call setup_workspace_dogleg(n,m,workspace%calculate_step_ws%dogleg_ws, & 
               options, inform)
          if (inform%alloc_status > 0) goto 1010

       case(2) ! use the AINT method
          call setup_workspace_AINT_tr(n,m,workspace%calculate_step_ws%AINT_tr_ws, & 
               options, inform)
          if (inform%alloc_status > 0) goto 1010

       case(3) ! More-Sorensen 
          call setup_workspace_more_sorensen(n,m,&
               workspace%calculate_step_ws%more_sorensen_ws,options,inform)
          if (inform%alloc_status > 0) goto 1010

       case (4) ! dtrs (Galahad)
          call setup_workspace_solve_galahad(n,m, & 
               workspace%calculate_step_ws%solve_galahad_ws, options, inform)
          if (inform%alloc_status > 0) goto 1010
       end select

       ! evaluate model in the main routine...       
       call setup_workspace_evaluate_model(n,m,workspace%evaluate_model_ws,options,inform)
       if (inform%alloc_status > 0) goto 1010

       return

       ! Error statements
1000   continue ! bad allocation from this subroutine
       inform%bad_alloc = 'setup_workspaces'
       inform%status = ERROR%ALLOCATION
       return

1010   continue ! bad allocation from called subroutine
       return

     end subroutine setup_workspaces

     subroutine remove_workspaces(workspace,options)

       type( NLLS_workspace ), intent(out) :: workspace
       type( nlls_options ), intent(in) :: options

       workspace%first_call = 0

       if(allocated(workspace%y)) deallocate(workspace%y)
       if(allocated(workspace%y_sharp)) deallocate(workspace%y_sharp)
       if(allocated(workspace%g_old)) deallocate(workspace%g_old)
       if(allocated(workspace%g_mixed)) deallocate(workspace%g_mixed)

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

       select case (options%nlls_method)

       case (1) ! use the dogleg method
          call remove_workspace_dogleg(workspace%calculate_step_ws%dogleg_ws, & 
               options)

       case(2) ! use the AINT method
          call remove_workspace_AINT_tr(workspace%calculate_step_ws%AINT_tr_ws, & 
               options)

       case(3) ! More-Sorensen 
          call remove_workspace_more_sorensen(&
               workspace%calculate_step_ws%more_sorensen_ws,options)

       case (4) ! dtrs (Galahad)
          call remove_workspace_solve_galahad(& 
               workspace%calculate_step_ws%solve_galahad_ws, options)

       end select

       ! evaluate model in the main routine...       
       call remove_workspace_evaluate_model(workspace%evaluate_model_ws,options)

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
       deallocate(w%work)
       allocate(w%work(lwork))     

       return

       ! Error statements
9000   continue  ! bad allocations in this subroutine
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

       return

     end subroutine remove_workspace_get_svd_J

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

       return

       ! Error statements
9000   continue  ! bad allocations in this subroutine
       inform%bad_alloc = 'setup_workspace_dogleg'
       inform%status = ERROR%ALLOCATION
       return

9010   continue  ! bad allocations from dependent subroutine
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

       return

9000   continue  ! local allocation error
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

       return

     end subroutine remove_workspace_solve_LLS

     subroutine setup_workspace_evaluate_model(n,m,w,options,inform)
       integer, intent(in) :: n, m        
       type( evaluate_model_work ) :: w
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(out) :: inform

       allocate( w%Jd(m), stat = inform%alloc_status )
       if (inform%alloc_status > 0) goto 9000
       allocate( w%Hd(n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000

       return

9000   continue
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = 'evaluate_model'
       return
     end subroutine setup_workspace_evaluate_model

     subroutine remove_workspace_evaluate_model(w,options)
       type( evaluate_model_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%Jd )) deallocate( w%Jd ) 
       if(allocated( w%Hd )) deallocate( w%Hd ) 

       return

     end subroutine remove_workspace_evaluate_model

     subroutine setup_workspace_AINT_tr(n,m,w,options,inform)
       integer, intent(in) :: n, m 
       type( AINT_tr_work ) :: w
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(out) :: inform

       allocate(w%A(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%v(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
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

       return

9000   continue ! local allocation error
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "AINT_tr"
       !call allocation_error(options,'AINT_tr')
       return

9010   continue ! allocation error from called subroutine
       return

     end subroutine setup_workspace_AINT_tr

     subroutine remove_workspace_AINT_tr(w,options)
       type( AINT_tr_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%A )) deallocate(w%A)
       if(allocated( w%v )) deallocate(w%v)
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
       deallocate(workquery)
       allocate( w%work(lwork), stat = inform%alloc_status )
       if (inform%alloc_status > 0) goto 9000

       return

9000   continue      
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "min_eig_symm"
       return

9010   continue
       inform%status = ERROR%FROM_EXTERNAL
       inform%external_name = "lapack_dsyev"

9020   continue
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
       deallocate(workquery)
       allocate( w%work(lwork), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate( w%nullindex(2*n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate( w%vecisreal(2*n), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000

       return

9000   continue
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "max_eig"
       return

9020   continue
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

       return

9000   continue ! allocation error
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "solve_general"
       return

     end subroutine setup_workspace_solve_general

     subroutine remove_workspace_solve_general(w, options)
       type( solve_general_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%A )) deallocate( w%A ) 
       if(allocated( w%ipiv )) deallocate( w%ipiv ) 
       return

     end subroutine remove_workspace_solve_general

     subroutine setup_workspace_solve_galahad(n,m,w,options,inform)
       integer, intent(in) :: n,m
       type( solve_galahad_work ) :: w
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(out) :: inform

       allocate(w%A(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%ev(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%v(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%v_trans(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%ew(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%d_trans(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000

       call setup_workspace_all_eig_symm(n,m,w%all_eig_symm_ws,options,inform)
       if (inform%status > 0) goto 9010

       if (options%scale > 0) then
          call setup_workspace_apply_scaling(n,m,w%apply_scaling_ws,options,inform)
          if (inform%status > 0) goto 9010
       end if

       return

9000   continue ! allocation error here
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "solve_galahad"
       return

9010   continue  ! allocation error from called subroutine
       return

     end subroutine setup_workspace_solve_galahad

     subroutine remove_workspace_solve_galahad(w,options)
       type( solve_galahad_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%A )) deallocate(w%A)
       if(allocated( w%ev )) deallocate(w%ev)
       if(allocated( w%v )) deallocate(w%v)
       if(allocated( w%v_trans )) deallocate(w%v_trans)
       if(allocated( w%ew )) deallocate(w%ew)
       if(allocated( w%d_trans )) deallocate(w%d_trans)

       call remove_workspace_all_eig_symm(w%all_eig_symm_ws,options)
       if (options%scale > 0) then
          call remove_workspace_apply_scaling(w%apply_scaling_ws,options)
       end if

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
       deallocate(workquery)
       allocate( w%work(lwork), stat = inform%alloc_status )
       if (inform%alloc_status > 0) goto 8000

       return

8000   continue  ! allocation error
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "all_eig_sym"
       return

9000   continue ! error from lapack
       inform%status = ERROR%FROM_EXTERNAL
       inform%external_name = "lapack_dsyev"
       return

     end subroutine setup_workspace_all_eig_symm

     subroutine remove_workspace_all_eig_symm(w,options)
       type( all_eig_symm_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%work )) deallocate( w%work ) 

     end subroutine remove_workspace_all_eig_symm

     subroutine setup_workspace_more_sorensen(n,m,w,options,inform)
       integer, intent(in) :: n,m
       type( more_sorensen_work ) :: w
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(out) :: inform

       allocate(w%A(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%LtL(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%v(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%q(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%y1(n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000
       allocate(w%AplusSigma(n,n),stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 9000

       call setup_workspace_min_eig_symm(n,m,w%min_eig_symm_ws,options,inform)
       if (inform%status > 0) goto 9010

       if (options%scale > 0) then
          call setup_workspace_apply_scaling(n,m,w%apply_scaling_ws,options,inform)
          if (inform%status > 0) goto 9010
       end if

       return

9000   continue ! allocation error here
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "more_sorenesen"
       return

9010   continue ! error from called subroutine
       return

     end subroutine setup_workspace_more_sorensen

     subroutine remove_workspace_more_sorensen(w,options)
       type( more_sorensen_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%A )) deallocate(w%A)
       if(allocated( w%LtL )) deallocate(w%LtL)
       if(allocated( w%v )) deallocate(w%v)
       if(allocated( w%q )) deallocate(w%q)
       if(allocated( w%y1 )) deallocate(w%y1)
       if(allocated( w%AplusSigma )) deallocate(w%AplusSigma)

       call remove_workspace_min_eig_symm(w%min_eig_symm_ws,options)
       if (options%scale > 0) then
          call remove_workspace_apply_scaling(w%apply_scaling_ws,options)
       end if

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

       return

1000   continue ! allocation error here
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "apply_scaling"
       return

1010   continue ! error from lower down subroutine
       return

     end subroutine setup_workspace_apply_scaling

     subroutine remove_workspace_apply_scaling(w,options)
       type( apply_scaling_work ) :: w
       type( nlls_options ), intent(in) :: options

       if(allocated( w%diag )) deallocate( w%diag )
       if(allocated( w%ev )) deallocate( w%ev ) 

       return 

     end subroutine remove_workspace_apply_scaling



   end module ral_nlls_internal

