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
    integer :: n
    integer :: m
    logical :: allocated
  end type matrix

  type, extends(matrix) :: dense_matrix
    real(wp), allocatable :: data(:,:)
  end type dense_matrix

  private
  public :: matrix, dense_matrix, init_dense_matrix, free_dense_matrix, mult_mv 

contains

  subroutine init_dense_matrix(A, n, m, inform)
    ! Create a dense_matrix structure, optionally with initial data
    type(dense_matrix), intent(out) :: A
    integer, intent(in) :: n, m

    A%n = n
    A%m = m
    if (A%allocated == .false.) then
      allocate(A%data(m,n), stat=inform%alloc_status)
      A%allocated = .true.
      if (present(init_data)) then
        call PREC(lacpy)('A', m, n, init_data, m, A%data, m)
      end if 
    else
      write(*,*) "Warning: init_dense_matrix called on already allocated matrix."
    end if
  end subroutine init_dense_matrix

  subroutine free_dense_matrix(A)
    ! Free the data of a dense_matrix structure
    type(dense_matrix), intent(inout) :: A

    if (A%allocated == .true.) deallocate(A%data)
    A%allocated = .false.
  end subroutine free_dense_matrix

  subroutine init_dense_nocopy(input_data, A, n, m)
    ! Initializes a dense_matrix structure by 
    ! taking ownership of some input data without copying
    real(wp), intent(in) :: input_data(:,:)
    type(dense_matrix), intent(out) :: A
    integer, intent(in) :: n, m

    A%n = n
    A%m = m
    call move_alloc(input_data, A%data)
  end subroutine init_dense_nocopy

  ! compute y = alpha * A * x + beta * y, where A is a dense matrix
  ! `trans` specifies whether to use A or its transpose
  ! this is an abstraction of the BLAS routine `gemv`
  subroutine mult_mv_dense(trans, A, x, y, alpha, beta)
    character(len=1), intent(in) :: trans
    type(dense_matrix), intent(in) :: A
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: alpha, beta

    call PREC(gemv)(trans, m, n, alpha, A%data, A%nrows, x, 1, beta, y, 1)
    
  end subroutine mult_mv_dense

end module MODULE_PREC(ral_nlls_matrix)
