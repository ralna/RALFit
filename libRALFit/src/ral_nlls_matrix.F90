! Copyright (c) 2020, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_matrix :: matrix abstractions

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_matrix)

  use MODULE_PREC(ral_nlls_workspaces), only: wp

  implicit none

  interface mult_mv 
    module procedure mult_mv_dense 
  end interface mult_mv

  type :: matrix
    integer :: nrows
    integer :: ncols
  end type matrix

  type, extends(matrix) :: dense_matrix
    real(wp), pointer :: data(:,:) => null()
  end type dense_matrix

contains

  ! compute y = alpha * A * x + beta * y, where A is a dense matrix
  ! `trans` specifies whether to use A or its transpose
  ! this is an abstraction of the BLAS routine `gemv`
  subroutine mult_mv_dense(trans, A, x, y, alpha, beta)
    character(len=1), intent(in) :: trans
    type(dense_matrix), intent(in) :: A
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: alpha, beta

    call PREC(gemv)(trans, A%nrows, A%ncols, alpha, A%data, A%nrows, x, 1, beta, y, 1)
    
  end subroutine mult_mv_dense

end module MODULE_PREC(ral_nlls_matrix)
