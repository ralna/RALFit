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
                                              NLLS_ERROR_BAD_SKETCH_METHOD, NLLS_ERROR_SPARSE_UNSUPPORTED
  use MODULE_PREC(lsqr_reverse), only: lsqr, lsqr_options, lsqr_inform, lsqr_keep, lsqr_free
  use MODULE_PREC(dct_module), only: dct, dct1d 
  use MODULE_PREC(ral_nlls_matrix), only: matrix, dense_matrix

  implicit none

  private
  public :: solve_LLS

  ! this interface should be REMOVED once all algorithms use the matrix interface
  interface solve_LLS
    module procedure solve_LLS
    module procedure solve_LLS_workaround
  end interface

  interface solve_direct
    module procedure solve_direct_dense
  end interface solve_direct

  ! note it is important that implementations of solve_cholesky
  ! throw an error when A is not PD! algorithms rely on this 
  ! for testing positive-definiteness of A
  interface solve_cholesky
    module procedure solve_cholesky_dense
  end interface solve_cholesky

contains
  ! WORKAROUND for while not all code is using the matrix interface
  subroutine solve_LLS_workaround(A, b, n, m, inform, w, options, pd)
    implicit none
    real(wp), intent(inout), contiguous :: A(:,:), b(:)
    integer, intent(in) :: n, m
    type(NLLS_inform), intent(inout) :: inform
    type(solve_LLS_work), intent(inout) :: w
    type(NLLS_options), intent(in) :: options 
    logical, intent(in) :: pd

    type(dense_matrix) :: A_mat

    call A_mat%init_matrix(n, m, options%fortran_jacobian, inform%alloc_status)
    call PREC(lacpy)('A', m, n, A, m, A_mat%data, m)
    call solve_LLS(A_mat, b, n, m, inform, w, options, pd)
    call PREC(lacpy)('A', m, n, A_mat%data, m, A, m)
    call A_mat%free_matrix()

  end subroutine solve_LLS_workaround


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
    class( matrix ), intent(inout) :: A
    real(wp), intent(inout), contiguous :: b(:)
    integer, intent(in) :: n, m
    type(NLLS_inform), intent(inout) :: inform
    type(solve_LLS_work), intent(inout) :: w
    type(NLLS_options), intent(in) :: options 
    logical, intent(in) :: pd

    ! if A is positive definite, we need to use Cholesky solver
    ! as routines rely on its test for positive-definiteness
    if (pd) then
      select type (A)
      type is (dense_matrix)
        call solve_cholesky(A, b, n, inform)
      class default
        inform%status = NLLS_ERROR_SPARSE_UNSUPPORTED
      end select
      return
    end if

    select case (options%lls_solver)
      case (1)  ! LAPACK
        select type (A)
        type is (dense_matrix)
          call solve_direct(A, b, n, m, inform, w, options)
        class default
          inform%status = NLLS_ERROR_SPARSE_UNSUPPORTED
        end select
      case (2)  ! use LSQR
        if (.not. w%lsqr_ws%allocated) then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          return
        end if
        call solve_lsqr(A, b, n, m, inform, w%lsqr_ws, options)
      case (3)  ! Randomised method (sketch-and-precondition)
        if (.not. w%lsqr_ws%allocated .or. .not. w%rand_ws%allocated) then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          return
        end if
        call solve_rand(A, b, n, m, options%sketch_size, inform, w, options)
      case default
        inform%status = NLLS_ERROR_BAD_LLS_SOLVER
    end select
    ! if LSQR or randomised fails, fallback to LAPACK GELS
    if (inform%status /= 0 .and. options%lls_solver > 1 .and. options%allow_fallback_method) then
      inform%status = 0
      select type (A)
        type is (dense_matrix)
          call solve_direct(A, b, n, m, inform, w, options)
        class default
          inform%status = NLLS_ERROR_SPARSE_UNSUPPORTED
      end select
    end if

  end subroutine solve_LLS

  subroutine solve_direct_dense(A, b, n, m, inform, w, options)
  ! Solve the least squares problem directly using LAPACK.
  type( dense_matrix ), intent(inout) :: A
  real(wp), intent(inout), contiguous :: b(:)
  integer, intent(in) :: n, m
  type(NLLS_inform), intent(inout) :: inform
  type(solve_LLS_work), intent(inout) :: w
  type(NLLS_options), intent(in) :: options

  integer :: lwork, lda, ldb

  if (.not. w%allocated) then
    inform%status = NLLS_ERROR_WORKSPACE_ERROR
    return
  end if

  if (n == m) then
    call PREC(gesv)(n, 1, A%data, n, w%ipiv, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?gesv'
    end if

  else
    lwork = size(w%work)
    if (options%Fortran_Jacobian) then
      lda = m
      ldb = max(m,n)
      call PREC(gels)('N', m, n, 1, A%data, lda, b, ldb, w%work, lwork, &
           inform%external_return)
    else
      lda = n
      ldb = max(m,n)
      call PREC(gels)('T', n, m, 1, A%data, lda, b, ldb, w%work, lwork, &
           inform%external_return)
    end if
  end if
  if (inform%external_return /= 0 ) then
    inform%status = NLLS_ERROR_FROM_EXTERNAL
    inform%external_name = 'lapack_?gels'
  end if
  end subroutine solve_direct_dense

  subroutine solve_cholesky_dense(A,b,n,inform)
    ! Wrapper around LAPACK's ?posv for positive definite systems
    implicit none
    type ( dense_matrix ), intent(inout) :: A
    REAL(wp), DIMENSION(:), INTENT(INOUT), contiguous :: b
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform

    call PREC(posv)('L', n, 1, A%data, n, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?posv'
    end if

  end subroutine solve_cholesky_dense 

  subroutine solve_lsqr(A, b, n, m, inform, w, options, precon, init_guess)
  ! non-preconditioned LSQR iterative solver for large sparse systems
  ! Implementation by Jennifer Scott, using the Papez-Tichy stopping criterion
  ! This is a reverse communication implementation, so we need
  ! a wrapper to handle the matrix products ourselves 
  ! Optional arguments:
  !   precon: if present, a preconditioner to apply in the matrix products
  !   init_guess: if true, assume w%x has been initialised with an initial guess
  ! currently, we assume the preconditioner is upper-triangular of size n x n
  integer, intent(in) :: m, n
  class( matrix ), intent(inout) :: A
  real(wp), intent(inout), contiguous :: b(:)
  type(NLLS_inform), intent(inout) :: inform
  type(LLS_lsqr_work), intent(inout) :: w
  type(NLLS_options), intent(in) :: options
  real(wp), intent(in), optional, contiguous :: precon(:,:)
  logical, intent(in), optional :: init_guess

  integer :: action, use_precon, rows, cols
  character(len=1) :: trans, no_trans
  real(wp) :: norm_a
  type(lsqr_options) :: lsqr_opt

  lsqr_opt%stop_test = 2
  lsqr_opt%itnlim = 5 * n

  ! stopping criterion requires an estimate of the spectral norm of A, 
  ! which we can get from a few iterations of the power method
  call estimate_norm(A, norm_a, 15, w, options%fortran_jacobian)

  if (options%fortran_jacobian) then
    rows = m
    cols = n
    trans = 'T'
    no_trans = 'N'
  else
    rows = n
    cols = m
    trans = 'N'
    no_trans = 'T'
  end if
  
  w%u = b
  w%v = 0.0_wp
  w%z = 0.0_wp
  w%ATu = 0.0_wp

  if (present(init_guess)) then
    if (init_guess) then
      ! update `u` for LSQR to be b - A x_0, so that we can start LSQR with this initial guess
      call A%mult_mv('N', w%x, w%u, -1.0_wp, 1.0_wp)
    else
      w%x = 0.0_wp
    end if
  else
    w%x = 0.0_wp
  end if

  action = 0
  if (present(precon)) then
    use_precon = 1
  else
    use_precon = 0
  end if

  ! first call to lsqr initializes the algorithm and returns the first action
  call lsqr(use_precon, action, m, n, w%u, w%v, w%z, w%x, w%keep, lsqr_opt, w%inform, anorm=norm_a)
 
  if (use_precon == 0) then
    do while (action /= 0)
      select case (action)
        case (1)  ! compute v = v + A^T u
          call A%mult_mv('T', w%u, w%v, 1.0_wp, 1.0_wp)
        case (2)  ! compute u = u + A v
          call A%mult_mv('N', w%v, w%u, 1.0_wp, 1.0_wp)
        end select
      call lsqr(0, action, m, n, w%u, w%v, w%z, w%x, w%keep, lsqr_opt, w%inform, anorm=norm_a)
    end do
  else
    do while (action /= 0)
      select case (action)
        case (1)  ! compute v = v + P^{-1} A^T u
          ! we don't want to explicitly invert R (for numerical stability), so we do this in two steps:
          ! first compute ATu = A^T u, then let t = R^{-1} ATu (i.e. solve R t = ATu)
          call A%mult_mv('T', w%u, w%ATu, 1.0_wp, 0.0_wp)
          call PREC(trtrs)('U', 'N', 'N', n, 1, precon, n, w%ATu, n, inform%external_return)
          if (inform%external_return /= 0) then
            inform%status = NLLS_ERROR_FROM_EXTERNAL
            inform%external_name = 'lapack_?trtrs'
            return
          end if
          w%v = w%v + w%ATu 
        case (2)  ! compute u = u + Az
          call A%mult_mv('N', w%z, w%u, 1.0_wp, 1.0_wp)
        case (3)  ! compute z = P^{-T} v (or equivalently, solve R^T z = v)
          w%z = w%v
          call PREC(trtrs)('U', 'T', 'N', n, 1, precon, n, w%z, n, inform%external_return)
          if (inform%external_return /= 0) then
            inform%status = NLLS_ERROR_FROM_EXTERNAL
            inform%external_name = 'lapack_?trtrs'
            return
          end if
      end select
      call lsqr(1, action, m, n, w%u, w%v, w%z, w%x, w%keep, lsqr_opt, w%inform, anorm=norm_a)
    end do
  end if

  if (w%inform%flag /= 0) then
    inform%status = NLLS_ERROR_FROM_EXTERNAL
    inform%external_return = w%inform%flag
    inform%external_name = 'lsqr'
    return
  end if


  b(1:n) = w%x(1:n)

  end subroutine solve_lsqr

  subroutine solve_rand(A, b, n, m, sketch_size, inform, w, options)
  ! Sketch-and-precondition randomised solver with sketch-and-solve initialisation
  ! Sources:
  !   P. Avron, P. Maymounkov, S. Toledo, 
  !   "Blendenpik: Supercharging LAPACK's Least-Squares Solver"
  !   URL: https://pdos.csail.mit.edu/~petar/papers/blendenpik-v1.pdf
  !   Meier, M., Nakatsukasa, Y., Townsend, A., Webb, M.
  !   "Are sketch-and-precondition linear solvers numerically stable?"
  !   URL: https://arxiv.org/abs/2302.07202
    class( matrix ), intent(inout) :: A
    real(wp), intent(inout), contiguous :: b(:)
    integer, intent(in) :: sketch_size, m, n
    type(NLLS_inform), INTENT(INOUT) :: inform
    type(solve_LLS_work), Intent(inout) :: w
    type(NLLS_options), Intent(In) :: options
    type(lsqr_options) :: lsqr_opt

    integer :: i, lwork, action, tsize, s
    real(wp) :: norm_a

    lwork = size(w%work)
    lsqr_opt%stop_test = 2

    w%rand_ws%temp_1 = 0.0_wp
    w%rand_ws%R = 0.0_wp
    w%rand_ws%SM = 0.0_wp
    w%rand_ws%Sb = 0.0_wp

    if (sketch_size == -1) then  ! default
      s = 4 * n
    else if (sketch_size <= n .or. sketch_size > m) then
      write(*,*) "Sketch size must be between n and m; setting default of 4n"
      s = 4 * n
    else
      s = sketch_size
    end if

    ! get sketch matrix SM and Sb
    call sketch(A, b, m, n, s, w%rand_ws%SM, w%rand_ws%Sb, &
                options, w%rand_ws, inform)
    if (inform%status /= 0) then
      return
    end if

    ! now take QR decomposition of SM
    call PREC(geqrf)(s, n, w%rand_ws%SM, s, w%rand_ws%temp_1, w%work, lwork, inform%external_return)
    if (inform%external_return /= 0) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?geqrf'
      return
    end if

    ! extract R from the output of geqrf (the upper triangle of the output matrix)
    do i = 1, n
      w%rand_ws%R(i, i:n) = w%rand_ws%SM(i, i:n)
    end do

    ! calculate initial guess from sketch-and-solve: x_0 = R^{-1} Q^T Sb
    call PREC(ormqr)('L', 'T', s, 1, n, w%rand_ws%SM, s, w%rand_ws%temp_1, &
              w%rand_ws%Sb, s, w%work, lwork, inform%external_return)
    if (inform%external_return /= 0) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?ormqr'
      return
    end if
    call PREC(trtrs)('U', 'N', 'N', n, 1, w%rand_ws%R, n, w%rand_ws%Sb, n, inform%external_return)
    if (inform%external_return /= 0) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?trtrs'
      return
    end if
    w%lsqr_ws%x(1:n) = w%rand_ws%Sb(1:n)

    ! now use this initial guess in LSQR to solve the original system, using R as a preconditioner 
    call solve_lsqr(A, b, n, m, inform, w%lsqr_ws, options, precon=w%rand_ws%R, init_guess=.true.)

  end subroutine solve_rand 

  subroutine estimate_norm(A, norm_a, num_iters, w, fortran_jacobian)
  ! power method to estimate spectral norm of A
    class(matrix), intent(in) :: A
    integer, intent(in) :: num_iters
    real(wp), intent(out) :: norm_a
    type(LLS_lsqr_work), intent(inout) :: w
    logical, intent(in) :: fortran_jacobian

    integer :: i

    ! start with a random vector
    call random_number(w%v)

    ! each iteration of the power method, we multiply by A^T A and renormalize
    ! this is done in two steps to be better-conditioned
    do i = 1, num_iters
      w%v = w%v/norm2(w%v)
      call A%mult_mv('N', w%v, w%u, 1.0_wp, 0.0_wp)
      call A%mult_mv('T', w%u, w%v, 1.0_wp, 0.0_wp)
    end do

    norm_a = sqrt(norm2(w%v))

  end subroutine estimate_norm

  subroutine sketch(A, b, m, n, sketch_size, SA, Sb, options, w, inform)
    ! Given a matrix A, return a "sketch" SA of size `sketch_size x n`,
    ! as well as the corresponding sketch Sb of b.
    class( matrix ), intent(in) :: A
    real(wp), intent(inout), contiguous :: b(:)
    integer, intent(in) :: m, n, sketch_size
    real(wp), intent(out), contiguous :: SA(:,:), Sb(:)
    type(NLLS_options), intent(in) :: options
    type(LLS_rand_work), intent(inout) :: w
    type(NLLS_inform), intent(inout) :: inform

    select case (options%sketch_method)
      case (1)  ! random sampling of rows
        call select_random_rows(A, b, n, sketch_size, SA, Sb, w)
      case (2)  ! random projection using DCT
        select type (A)
        type is (dense_matrix)
          call dct_sketch(A, b, m, n, sketch_size, SA, Sb, w)
        class default
          inform%status = NLLS_ERROR_SPARSE_UNSUPPORTED
        end select
      case default
        inform%status = NLLS_ERROR_BAD_SKETCH_METHOD
    end select

  end subroutine sketch

  subroutine select_random_rows(A, b, n, s, SA, Sb, w)
    ! Select `s` random rows of A to form SA
    class( matrix ), intent(in) :: A
    real(wp), intent(inout), contiguous :: b(:)
    integer, intent(in) :: n, s
    real(wp), intent(out) :: SA(:,:), Sb(:)
    type(LLS_rand_work), intent(inout) :: w

    integer :: i, idx

    call random_number(w%temp_2)

    ! to select `s` random rows, we generate `m` random numbers in [0,1),
    ! and take the indices of the smallest `s` numbers
    do i = 1, s
      idx = minloc(w%temp_2, dim=1)
      w%temp_2(idx) = 1.0_wp
      call A%copy_row(1.0_wp, idx, SA(i, :)) 
      Sb(i) = b(idx)
    end do

  end subroutine select_random_rows

  subroutine dct_sketch(A, b, m, n, s, SA, Sb, w)
    ! Create a sketch of A by taking a random subsample of rows
    ! applying random sign flips and a DCT to reduce coherence
    ! (i.e. the probability that the sketch is rank deficient)
    class( dense_matrix ), intent(in) :: A
    real(wp), intent(inout), contiguous :: b(:)
    integer, intent(in) :: m, n, s 
    real(wp), intent(inout) :: SA(:,:), Sb(:)
    type(LLS_rand_work), intent(inout) :: w

    integer :: i 
    real(wp) :: rand_no

    ! TODO: generalise this to a general matrix rather than a dense one

    w%DCT_A%data = 0.0_wp

    ! create DA, which is A with random sign flips
    ! i.e. D is a random diagonal matrix with +/- 1 on the diagonal
    do i = 1, m
      call random_number(rand_no)
      if (rand_no < 0.5_wp) then
        call A%copy_row(-1.0_wp, i, w%DCT_A%data(i, :))
        w%temp_1(i) = -b(i)
      else
        call A%copy_row(1.0_wp, i, w%DCT_A%data(i, :))
        w%temp_1(i) = b(i)
      end if
    end do

    ! apply discrete cosine transform to each column
    call dct(w%DCT_A%data, n, m, w%dct_work)
    call dct1d(w%temp_1, m, w%dct_work)

    ! and normalise DCT
    w%DCT_A%data = w%DCT_A%data / sqrt(2.0_wp * real(m, wp))
    w%temp_1 = w%temp_1 / sqrt(2.0_wp * real(m, wp))

    ! and finaly select `s` random rows of DCT_A to form SA
    ! and the corresponding entries of b to form Sb 
    call select_random_rows(w%DCT_A, w%temp_1, n, s, SA, Sb, w)

  end subroutine dct_sketch

end module MODULE_PREC(ral_nlls_linear)
