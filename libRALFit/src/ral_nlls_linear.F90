! Copyright (c) 2020, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_linear :: linear solvers for internal use in RALFit

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_linear)
  use MODULE_PREC(ral_nlls_workspaces), only: wp, solve_LLS_work, NLLS_inform, &
                                              NLLS_options, NLLS_ERROR_WORKSPACE_ERROR, &
                                              NLLS_ERROR_FROM_EXTERNAL, NLLS_ERROR_BAD_LLS_SOLVER
  implicit none

  private
  public :: solve_LLS, solve_LLS_nocopy

contains
  ! todo: Jacobian argument
  subroutine solve_LLS(A,b,A_out,x,n,m,inform,w,options,pd)
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------
!  Solves the linear least squares problem Ax = b
!  using the method specified in options%lls_solver
!  Input:
!    A         - LHS matrix of the least squares problem 
!    b         - The RHS, overwritten with x on output 
!    n         - Number of columns in A
!    m         - Number of rows in A
!    inform    - NLLS_inform structure
!    options   - NLLS_options structure
!    w         - workspace structure
!    pd        - logical, true if A is known to be positive definite
!  Output:
!    A_out     - output copy of the original A matrix
!    x         - solution vector 
!    inform    - NLLS_inform structure
!    work      - updated workspace structure
!  -----------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, m
    real(wp), intent(in) :: A(n,m)
    real(wp), intent(in) :: b(m)
    real(wp), intent(out) :: A_out(n,m)
    real(wp), intent(out) :: x(n)
    type(NLLS_inform), intent(inout) :: inform
    type(NLLS_options), Intent(In) :: options 
    type(solve_LLS_work), intent(inout) :: w
    logical, intent(in) :: pd

    real(wp) :: b_out(m)

    A_out = A
    b_out = b

    call solve_LLS_nocopy(A_out, b_out, n, m, inform, w, options, pd)
    if (inform%status == 0) then
      x(1:n) = b_out(1:n)
    end if

  end subroutine solve_LLS

  subroutine solve_LLS_nocopy(A, b, n, m, inform, w, options, pd)
!  linear solver core which overwrites A and b. 
    implicit none
    integer, intent(in) :: n, m
    real(wp), dimension(m,n), intent(inout) :: A
    real(wp), dimension(m), intent(inout) :: b
    type(NLLS_inform), intent(inout) :: inform
    type(NLLS_options), intent(in) :: options 
    type(solve_LLS_work), intent(inout) :: w
    logical, intent(in) :: pd

    select case (options%lls_solver) 
      case (1)  ! LAPACK
        if (pd) then
          call solve_posv(A,b,n,inform)
        else
          if (n == m) then
            call solve_gesv(A,b,n,inform)
          else
            if (.not. w%allocated) then
              inform%status = NLLS_ERROR_WORKSPACE_ERROR
              goto 100
            end if
            call solve_gels(A,b,n,m,inform,w,options)
          end if
        end if
      case default
        inform%status = NLLS_ERROR_BAD_LLS_SOLVER
    end select

100 continue
  end subroutine solve_LLS_nocopy

  subroutine solve_gels(A,b,n,m,inform,w,options)
!   Wrapper around LAPACK's ?gels
   implicit none 
   REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: A
   REAL(wp), DIMENSION(:), INTENT(INOUT) :: b
   INTEGER, INTENT(IN) :: n, m
   type(NLLS_inform), INTENT(INOUT) :: inform
   type(NLLS_options), Intent(In) :: options 

   integer :: lwork, lda, ldb
   type( solve_LLS_work ), Intent(inout) :: w

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

  subroutine solve_gesv(A,b,n,inform)
!   Wrapper around LAPACK's ?gesv
    implicit none
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: A
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform

    ! NB: we never actually use ipiv in the library
    integer, dimension(n) :: ipiv

    call PREC(gesv)(n, 1, A, n, ipiv, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?gesv'
    end if

  end subroutine solve_gesv

  subroutine solve_posv(A,b,n,inform)
!   Wrapper around LAPACK's ?posv for positive definite systems
    implicit none
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: A
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform

    call PREC(posv)('L', n, 1, A, n, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?posv'
    end if

  end subroutine solve_posv

end module MODULE_PREC(ral_nlls_linear)
