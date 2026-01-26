! ral_nlls_linear :: linear solvers for internal use in RALFit

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_linear)
  use MODULE_PREC(ral_nlls_workspaces), only: wp, solve_LLS_work, NLLS_inform, &
                                              NLLS_options, NLLS_ERROR_WORKSPACE_ERROR, &
                                              NLLS_ERROR_FROM_EXTERNAL
  implicit none

  private
  public :: solve_LLS, solve_LLS_nocopy

contains
  subroutine solve_LLS(A,b,A_out,x,n,m,inform,w,options,pd)
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------
!  Solves the linear least squares problem Ax = b
!  using the method specified in options%lls_solver
!  Input:
!    A         - LHS matrix of the least squares problem (m x n)
!    b         - The RHS, overwritten with x on output (m x 1)
!    n         - Number of columns in A
!    m         - Number of rows in A
!    inform    - NLLS_inform structure
!    options   - NLLS_options structure
!    w         - workspace structure
!    pd        - logical, true if A is known to be positive definite
!  Output:
!    A_out     - output copy of the original A matrix (m x n)
!    x         - solution vector (n x 1)
!    inform    - NLLS_inform structure
!    work      - updated workspace structure
!  -----------------------------------------------------------------
    implicit none
    real(wp), dimension(:), intent(in) :: A
    real(wp), dimension(:), intent(in) :: b
    real(wp), dimension(:), intent(out) :: A_out
    real(wp), dimension(:), intent(out) :: x
    integer, intent(in) :: n, m
    type(NLLS_inform), intent(inout) :: inform
    type(NLLS_options), Intent(In) :: options 
    type(solve_LLS_work), intent(inout) :: w
    logical, intent(in) :: pd

    A_out = A
    x = b

    call solve_LLS_nocopy(A_out, x, n, m, inform, w, options, pd)

  end subroutine solve_LLS

  subroutine solve_LLS_nocopy(A, b, n, m, inform, w, options, pd)
!  linear solver core which overwrites A and b. 
    implicit none
    real(wp), dimension(:), intent(inout) :: A
    real(wp), dimension(:), intent(inout) :: b
    integer, intent(in) :: n, m
    type(NLLS_inform), intent(inout) :: inform
    type(NLLS_options), intent(in) :: options 
    type(solve_LLS_work), intent(inout) :: w
    logical, intent(in) :: pd

    select case (options%lls_solver) 
      case (1)  ! LAPACK
        if (pd) then
          call solve_posv(A,b,n,inform,options)
        else
          if (.not. w%allocated) then
            inform%status = NLLS_ERROR_WORKSPACE_ERROR
            goto 100
          end if
          call solve_gels(A,b,n,m,inform,w,options)
        end if
    end select

    if (inform%external_return /= 0 ) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
    end if 

100   continue
  end subroutine solve_LLS_nocopy

  subroutine solve_gels(A,b,n,m,inform,w,options)
!   Wrapper around LAPACK's ?gels
   implicit none 
   REAL(wp), DIMENSION(:), INTENT(IN) :: A
   REAL(wp), DIMENSION(:), INTENT(IN) :: b
   INTEGER, INTENT(IN) :: n, m
   type(NLLS_inform), INTENT(INOUT) :: inform
   type(NLLS_options), Intent(In) :: options 

   integer, Parameter :: nrhs = 1
   integer :: lwork, lda, ldb
   type( solve_LLS_work ), Intent(inout) :: w

   lwork = size(w%work)

   if (options%Fortran_Jacobian) then
      lda = m
      ldb = max(m,n)
      call PREC(gels)('N', m, n, nrhs, A, lda, w%temp, ldb, w%work, lwork, &
           inform%external_return)
   else
      lda = n
      ldb = max(m,n)
      call PREC(gels)('T', n, m, nrhs, A, lda, w%temp, ldb, w%work, lwork, &
           inform%external_return)
   end if
   if (inform%external_return /= 0 ) then
      inform%status = NLLS_ERROR_FROM_EXTERNAL
      inform%external_name = 'lapack_?gels'
   end if

  end subroutine solve_gels

  subroutine solve_posv(A,b,n,inform,options)
!   Wrapper around LAPACK's ?posv for positive definite systems
    implicit none
    REAL(wp), DIMENSION(:), INTENT(IN) :: A
    REAL(wp), DIMENSION(:), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: n
    type(NLLS_inform), INTENT(INOUT) :: inform
    type(NLLS_options), Intent(In) :: options

    call PREC(posv)('N', n, 1, A, n, b, n, inform%external_return)
    if (inform%external_return /= 0 ) then
       inform%status = NLLS_ERROR_FROM_EXTERNAL
       inform%external_name = 'lapack_?posv'
    end if

  end subroutine solve_posv

end module MODULE_PREC(ral_nlls_linear)
