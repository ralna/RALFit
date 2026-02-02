! Copyright (c) 2020, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_linear :: linear solvers for internal use in RALFit

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_linear)
  use MODULE_PREC(ral_nlls_workspaces), only: wp, solve_LLS_work, LLS_lsqr_work, LLS_rand_work, &
                                              NLLS_inform, NLLS_options, NLLS_ERROR_WORKSPACE_ERROR, &
                                              NLLS_ERROR_FROM_EXTERNAL, NLLS_ERROR_BAD_LLS_SOLVER, &
                                              NLLS_ERROR_BAD_SKETCH_METHOD, NLLS_ERROR_BAD_SKETCH_SIZE
  use MODULE_PREC(lsqr_reverse), only: lsqr, lsqr_options, lsqr_inform, lsqr_keep, lsqr_free
  use MODULE_PREC(dct_module), only: dct 

  implicit none

  private
  public :: solve_LLS

contains
  ! todo: Jacobian argument
  subroutine solve_LLS(A, b, n, m, inform, w, options, pd)
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------
!  Solves the linear least squares problem Ax = b
!  using the method specified in options%lls_solver
!  Input:
!    A         - LHS matrix of the least squares problem 
!    b         - The RHS, overwritten with result x on output
!    n         - Number of columns in A
!    m         - Number of rows in A
!    inform    - NLLS_inform structure
!    w         - workspace structure
!    options   - NLLS_options structure
!    pd        - logical, true if A is known to be positive definite
!  -----------------------------------------------------------------
    implicit none
    real(wp), intent(inout), contiguous :: A(:,:), b(:)
    integer, intent(in) :: n, m
    type(NLLS_inform), intent(inout) :: inform
    type(solve_LLS_work), intent(inout) :: w
    type(NLLS_options), intent(in) :: options 
    logical, intent(in) :: pd

    select case (options%lls_solver)
      case (1)  ! LAPACK
        if (pd) then
          call solve_posv(A,b,n,inform)
        else
          if (n == m) then
            call solve_gesv(A,b,n,inform,w)
          else
            call solve_gels(A,b,n,m,inform,w,options)
          end if
        end if
      case (2)  ! use LSQR
        if (.not. w%lsqr_ws%allocated) then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          return
        end if
        call solve_lsqr(A, b, n, m, inform, w%lsqr_ws, options)
        ! with LAPACK fallback...
        if (inform%status /= 0) then
          if (options%allow_fallback_method) then
            inform%status = 0
            call solve_gels(A,b,n,m,inform,w,options)
          else
            return
          end if
        end if
      case (3)  ! Randomised method (Blendenpik)
        if (.not. w%lsqr_ws%allocated .or. .not. w%rand_ws%allocated) then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          return
        end if
        call blendenpik(A, b, n, m, options%sketch_size, inform, w, options)
        ! if Blendenpik fails, fall back to LAPACK
      case default
        inform%status = NLLS_ERROR_BAD_LLS_SOLVER
    end select
    ! if LSQR or Blendenpik fails, fallback to LAPACK GELS
    if (inform%status /= 0 .and. options%lls_solver > 1 .and. options%allow_fallback_method) then
      inform%status = 0
      call solve_gels(A,b,n,m,inform,w,options)
    end if

  end subroutine solve_LLS

  subroutine solve_gels(A,b,n,m,inform,w,options)
!   Wrapper around LAPACK's ?gels
   implicit none 
   real(wp), intent(inout), contiguous :: A(:,:), b(:)
   INTEGER, INTENT(IN) :: n, m
   type(NLLS_inform), INTENT(INOUT) :: inform
   type(NLLS_options), Intent(In) :: options 

   integer :: lwork, lda, ldb
   type( solve_LLS_work ), Intent(inout) :: w
   if (.not. w%allocated) then
     inform%status = NLLS_ERROR_WORKSPACE_ERROR
     return
   end if
   lwork = size(w%work)

   if (options%Fortran_Jacobian) then
      lda = m
      ldb = max(m,n)
      call PREC(gels)('N', m, n, 1, A, lda, b, ldb, w%work, lwork, &
           inform%external_return)
   else
      lda = n
      ldb = max(m,n)
      call PREC(gels)('T', n, m, 1, A, lda, b, ldb, w%work, lwork, &
           inform%external_return)
   end if
   if (inform%external_return /= 0 ) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?gels'
   end if

  end subroutine solve_gels

  subroutine solve_gesv(A,b,n,inform,w)
!   Wrapper around LAPACK's ?gesv
    implicit none
    real(wp), intent(inout), contiguous :: A(:,:), b(:)
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform
    type(solve_LLS_work), intent(inout) :: w

    call PREC(gesv)(n, 1, A, n, w%ipiv, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?gesv'
    end if

  end subroutine solve_gesv

  subroutine solve_posv(A,b,n,inform)
!   Wrapper around LAPACK's ?posv for positive definite systems
    implicit none
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), contiguous :: A
    REAL(wp), DIMENSION(:), INTENT(INOUT), contiguous :: b
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform

    call PREC(posv)('L', n, 1, A, n, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?posv'
    end if

  end subroutine solve_posv

  subroutine solve_lsqr(A, b, n, m, inform, w, options)
  ! non-preconditioned LSQR iterative solver for large sparse systems
  ! Implementation by Jennifer Scott, using the Papez-Tichy stopping criterion
  ! This is a reverse communication implementation, so we need
  ! a wrapper to handle the arithmetic ourselves
  integer, intent(in) :: m, n
  real(wp), intent(inout), contiguous :: A(:, :), b(:)
  type(NLLS_inform) :: inform
  type(LLS_lsqr_work), intent(inout) :: w
  type(NLLS_options), intent(in) :: options
  type(lsqr_options) :: lsqr_opt

  integer :: action
  real(wp) :: norm_a

  lsqr_opt%stop_test = 2
  lsqr_opt%itnlim = 5 * n

  ! stopping criterion requires an estimate of the spectral norm of A, 
  ! which we can get from a few iterations of the power method
  call estimate_norm(A, m, n, norm_a, 15, w)
  
  w%u = b
  w%v = 0.0_wp
  w%z = 0.0_wp
  w%x = 0.0_wp

  action = 0

  ! first call to lsqr initializes the algorithm and returns the first action
  call lsqr(0, action, m, n, w%u, w%v, w%z, w%x, w%keep, lsqr_opt, w%inform, anorm=norm_a)
  
  do while (action /= 0)
    select case (action)
      case (1)  ! compute v = v + A^T u
        call PREC(gemv)('T', m, n, 1.0_wp, A, m, w%u, 1, 1.0_wp, w%v, 1)
      case (2)  ! compute u = u + A v
        call PREC(gemv)('N', m, n, 1.0_wp, A, m, w%v, 1, 1.0_wp, w%u, 1)
      end select
    call lsqr(0, action, m, n, w%u, w%v, w%z, w%x, w%keep, lsqr_opt, w%inform, anorm=norm_a)
  end do

  if (w%inform%flag /= 0) then
    inform%status = NLLS_ERROR_FROM_EXTERNAL
    inform%external_return = w%inform%flag
    inform%external_name = 'lsqr'
    return
  end if


  b(1:n) = w%x(1:n)

  end subroutine solve_lsqr

  subroutine blendenpik(A, b, n, m, sketch_size, inform, w, options)
  ! Blendenpik (sketch-and-precondition) randomised solver
  ! Source:
  !   P. Avron, P. Maymounkov, S. Toledo, 
  !   "Blendenpik: Supercharging LAPACK's Least-Squares Solver"
  !   URL: https://pdos.csail.mit.edu/~petar/papers/blendenpik-v1.pdf
    real(wp), intent(inout), contiguous :: A(:,:), b(:)
    integer, intent(in) :: sketch_size, m, n
    type(NLLS_inform), INTENT(INOUT) :: inform
    type(solve_LLS_work), Intent(inout) :: w
    type(NLLS_options), Intent(In) :: options
    type(lsqr_options) :: lsqr_opt

    integer :: i, lwork, action
    real(wp) :: norm_a

    lwork = size(w%work)
    lsqr_opt%stop_test = 2

    w%rand_ws%tau = 0.0_wp
    w%rand_ws%R = 0.0_wp
    w%rand_ws%SM = 0.0_wp
    w%rand_ws%ATu = 0.0_wp
    w%rand_ws%DCT_A = 0.0_wp

    ! get sketch matrix SM
    call sketch(A, m, n, sketch_size, w%rand_ws%SM, options, w%rand_ws, inform)
    if (inform%status /= 0) then
      return
    end if

    ! now take QR decomposition of SM
    call PREC(geqrf)(sketch_size, n, w%rand_ws%SM, sketch_size, w%rand_ws%tau, w%work, lwork, inform%external_return)
    if (inform%external_return /= 0) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?geqrf'
      return
    end if

    ! todo: initial guess!
    ! extract R from QR decomposition
    do i = 1, n
      w%rand_ws%R(i, i:n) = w%rand_ws%SM(i, i:n)
    end do

    call estimate_norm(A, m, n, norm_a, 15, w%lsqr_ws)

    ! first, we run LSQR without preconditioning to 
    w%lsqr_ws%x = 0.0_wp

    ! solve using LSQR
    ! first call to lsqr initializes the algorithm and returns the first action
    w%lsqr_ws%u = b
    w%lsqr_ws%v = 0.0_wp
    w%lsqr_ws%z = 0.0_wp
    action = 0
    call lsqr(1, action, m, n, w%lsqr_ws%u, w%lsqr_ws%v, w%lsqr_ws%z, &
              w%lsqr_ws%x, w%lsqr_ws%keep, lsqr_opt, w%lsqr_ws%inform, &
              anorm=norm_a)
    do while (action /= 0)
      select case (action)
        case (1)  ! compute v = v + P^{-1} A^T u
          ! we don't want to explicitly invert R (for numerical stability), so we do this in two steps:
          ! first compute ATu = A^T u, then let t = R^{-1} ATu (i.e. solve R t = ATu)
          call PREC(gemv)('T', m, n, 1.0_wp, A, m, w%lsqr_ws%u, 1, 0.0_wp, w%rand_ws%ATu, 1)
          call PREC(trtrs)('U', 'N', 'N', n, 1, w%rand_ws%R, n, w%rand_ws%ATu, n, inform%external_return)
          if (inform%external_return /= 0) then
            inform%status = NLLS_ERROR_FROM_EXTERNAL
            inform%external_name = 'lapack_?trtrs'
            return
          end if
          w%lsqr_ws%v = w%lsqr_ws%v + w%rand_ws%ATu 
        case (2)  ! compute u = u + Az
          call PREC(gemv)('N', m, n, 1.0_wp, A, m, w%lsqr_ws%z, 1, 1.0_wp, w%lsqr_ws%u, 1)
        case (3)  ! compute z = P^{-T} v (or equivalently, solve R^T z = v)
          w%lsqr_ws%z = w%lsqr_ws%v
          call PREC(trtrs)('U', 'T', 'N', n, 1, w%rand_ws%R, n, w%lsqr_ws%z, n, inform%external_return)
          if (inform%external_return /= 0) then
            inform%status = NLLS_ERROR_FROM_EXTERNAL
            inform%external_name = 'lapack_?trtrs'
            return
          end if
      end select
      call lsqr(1, action, m, n, w%lsqr_ws%u, w%lsqr_ws%v, w%lsqr_ws%z, &
                w%lsqr_ws%x, w%lsqr_ws%keep, lsqr_opt, w%lsqr_ws%inform, &
                anorm=norm_a)
    end do
    if (w%lsqr_ws%inform%flag /= 0) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_return = w%lsqr_ws%inform%flag
      inform%external_name = 'lsqr'
      return
    end if

    b(1:n) = w%lsqr_ws%x(1:n) 

  end subroutine blendenpik

  subroutine estimate_norm(A, m, n, norm_a, num_iters, w)
  ! power method to estimate spectral norm of A
    real(wp), intent(in), contiguous :: A(:,:)
    integer, intent(in) :: m, n, num_iters
    real(wp), intent(out) :: norm_a
    type(LLS_lsqr_work), intent(inout) :: w

    integer :: i

    ! start with a random vector
    call random_number(w%v)

    ! each iteration of the power method, we multiply by A^T A and renormalize
    ! this is done in two steps to be better-conditioned
    do i = 1, num_iters
      w%v = w%v / norm2(w%v)
      call PREC(gemv)('N', m, n, 1.0_wp, A, m, w%v, 1, 0.0_wp, w%u, 1)
      call PREC(gemv)('T', m, n, 1.0_wp, A, m, w%u, 1, 0.0_wp, w%v, 1)
    end do

    norm_a = sqrt(norm2(w%v))

  end subroutine estimate_norm

  subroutine sketch(A, m, n, sketch_size, SA, options, w, inform)
    ! Given a matrix A, return a "sketch" SA of size `sketch_size x n`.
    real(wp), intent(in), contiguous :: A(:,:)
    integer, intent(in) :: m, n, sketch_size
    real(wp), intent(out), contiguous :: SA(:,:)
    type(NLLS_options), intent(in) :: options
    type(LLS_rand_work), intent(inout) :: w
    type(NLLS_inform), intent(inout) :: inform

    if (sketch_size <= 0 .or. sketch_size > m) then
      inform%status = NLLS_ERROR_BAD_SKETCH_SIZE
      return
    end if

    select case (options%sketch_method)
      case (1)  ! random sampling of rows
        call select_random_rows(A, n, sketch_size, SA, w)
      case (2)  ! random projection using DCT
        call dct_sketch(A, m, n, sketch_size, SA, w)
      case default
        inform%status = NLLS_ERROR_BAD_SKETCH_METHOD
    end select

  end subroutine sketch

  subroutine select_random_rows(A, n, s, SA, w)
    ! Select `s` random rows of A to form SA
    real(wp), intent(in), contiguous :: A(:,:)
    integer, intent(in) :: n, s
    real(wp), intent(out) :: SA(:,:)
    type(LLS_rand_work), intent(inout) :: w

    integer :: i, idx

    call random_number(w%temp)

    ! to select `s` random rows, we generate `m` random numbers in [0,1),
    ! and take the indices of the smallest `s` numbers
    do i = 1, s
      idx = minloc(w%temp, dim=1)
      w%temp(idx) = 1.0_wp
      SA(i, :) = A(idx, :)
    end do
  end subroutine select_random_rows

  subroutine dct_sketch(A, m, n, s, SA, w)
    ! Create a sketch of A by taking a random subsample of rows
    ! applying random sign flips and a DCT to reduce coherence
    ! (i.e. the probability that the sketch is rank deficient)
    real(wp), intent(in), contiguous :: A(:,:)
    integer, intent(in) :: m, n, s
    real(wp), intent(inout) :: SA(:,:)
    type(LLS_rand_work), intent(inout) :: w

    integer :: i
    real(wp) :: rand_no

    w%DCT_A = 0.0_wp

    ! create DA, which is A with random sign flips
    ! i.e. D is a random diagonal matrix with +/- 1 on the diagonal
    do i = 1, m
      call random_number(rand_no)
      if (rand_no < 0.5_wp) then
        w%DCT_A(i, :) = -A(i, :)
      else
        w%DCT_A(i, :) = A(i, :)
      end if
    end do

    ! apply discrete cosine transform to each column
    call dct(w%DCT_A, n, m)

    ! and normalise DCT
    w%DCT_A = w%DCT_A / sqrt(2.0_wp * real(m, wp))

    ! and finaly select `s` random rows of DCT_A to form SA
    call select_random_rows(w%DCT_A, n, s, SA, w)

  end subroutine dct_sketch

end module MODULE_PREC(ral_nlls_linear)
