! Copyright (c) 2020, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_matrix :: matrix abstractions

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_matrix)

  use MODULE_PREC(ral_nlls_types), only: wp

  implicit none

  type, abstract :: matrix
    integer :: n
    integer :: m
    logical :: allocated
    contains
      procedure(init_iface), deferred       :: init_matrix  ! initialise new matrix
      procedure(free_iface), deferred       :: free_matrix  ! free memory held by matrix
      procedure(mult_mv_iface), deferred    :: mult_mv      ! generalisation of *gemv
      procedure(copy_mat_iface), deferred   :: copy_matrix  ! generalisation of *lacpy
      procedure(copy_row_iface), deferred   :: copy_row     ! copy a row of a matrix to a vector
      procedure(mult_inner_iface), deferred :: mult_inner   ! compute A^T * A and store in out_A
  end type matrix

  abstract interface

    subroutine init_iface(A, n, m, fortran_order, status)
      import :: matrix 
      class(matrix), intent(inout) :: A
      integer, intent(in) :: n, m
      logical, intent(in) :: fortran_order
      integer, intent(out) :: status
    end subroutine init_iface

    subroutine free_iface(A, stat)
      import :: matrix 
      class(matrix), intent(inout) :: A
      integer, intent(out), optional :: stat
    end subroutine free_iface

    subroutine mult_mv_iface(A, trans, x, y, alpha, beta)
      import :: matrix, wp
      class(matrix), intent(in) :: A
      character(len=1), intent(in) :: trans
      real(wp), intent(in) :: x(:)
      real(wp), intent(inout) :: y(:)
      real(wp), intent(in) :: alpha, beta
    end subroutine mult_mv_iface

    subroutine copy_mat_iface(src, dest)
      import :: matrix, wp
      class(matrix), intent(in) :: src
      real(wp), intent(out) :: dest(:,:)
    end subroutine

    subroutine copy_row_iface(src, alpha, row, dest)
      import :: matrix, wp
      class(matrix), intent(in) :: src
      real(wp), intent(out) :: dest(:)
      integer, intent(in) :: row
      real(wp), intent(in) :: alpha
    end subroutine copy_row_iface

    subroutine mult_inner_iface(A, out_ATA)
      import :: matrix, wp
      class(matrix), intent(in) :: A
      real(wp), intent(out) :: out_ATA(:,:)
    end subroutine mult_inner_iface

  end interface

  type, extends(matrix) :: dense_matrix
    real(wp), allocatable :: data(:,:)
    logical :: fortran_order 
    contains
      procedure :: init_matrix => init_matrix_dense
      procedure :: free_matrix => free_matrix_dense
      procedure :: mult_mv => mult_mv_dense
      procedure :: copy_matrix => copy_matrix_dense
      procedure :: copy_row => copy_row_dense
      procedure :: mult_inner => mult_inner_dense
  end type dense_matrix

  private
  public :: matrix, dense_matrix

contains

  ! --------------------------------------------------------------------------------
  ! Dense matrix subroutines
  ! --------------------------------------------------------------------------------

  subroutine init_matrix_dense(A, n, m, fortran_order, status)
    ! Create a dense_matrix structure
    class(dense_matrix), intent(inout) :: A
    integer, intent(in) :: n, m
    logical, intent(in) :: fortran_order
    integer, intent(out) :: status

    A%n = n
    A%m = m
    A%fortran_order = fortran_order
    if (.not. allocated(A%data)) then
      if (fortran_order) then
        allocate(A%data(m, n))
      else
        allocate(A%data(n, m))
      end if
      A%allocated = .true.
    else
      write(*,*) "Warning: init_matrix_dense called on already allocated matrix."
    end if
  end subroutine init_matrix_dense 

  subroutine free_matrix_dense(A, stat)
    ! Free the data of a dense_matrix structure
    class(dense_matrix), intent(inout) :: A
    integer, intent(out), optional :: stat

    if (present(stat)) then
      if (A%allocated) deallocate(A%data, stat=stat)
    else
      if (A%allocated) deallocate(A%data)
    end if
    A%allocated = .false.
  end subroutine free_matrix_dense 

  ! compute y = alpha * A * x + beta * y, where A is a dense matrix
  ! `trans` specifies whether to use A or its transpose
  ! this is an abstraction of the BLAS routine `gemv`
  subroutine mult_mv_dense(A, trans, x, y, alpha, beta)
    class(dense_matrix), intent(in) :: A
    character(len=1), intent(in) :: trans
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: alpha, beta

    if (A%fortran_order) then
      call PREC(gemv)(trans, A%m, A%n, alpha, A%data, A%m, x, 1, beta, y, 1)
    else
      if (trans == 'T') then
        call PREC(gemv)('N', A%n, A%m, alpha, A%data, A%n, x, 1, beta, y, 1)
      else
        call PREC(gemv)('T', A%n, A%m, alpha, A%data, A%n, x, 1, beta, y, 1)
      end if
    end if
    
  end subroutine mult_mv_dense

  subroutine copy_matrix_dense(src, dest)
    ! Copy the contents of one dense_matrix to another
    class(dense_matrix), intent(in) :: src
    real(wp), intent(out) :: dest(:,:)

    call PREC(lacpy)('A', src%m, src%n, src%data, src%m, dest, src%m)
  end subroutine copy_matrix_dense

  subroutine copy_row_dense(src, alpha, row, dest)
    ! Copy a single row from one dense_matrix to a vector, applying a scalar 
    class(dense_matrix), intent(in) :: src
    real(wp), intent(out) :: dest(:)
    integer, intent(in) :: row
    real(wp), intent(in) :: alpha

    if (src%fortran_order) then
      dest = alpha * src%data(row, :)
    else
      dest = alpha * src%data(:, row)
    end if
  end subroutine copy_row_dense

  ! todo: we want A^T A to be a matrix type, but this will have a lot of knock-on
  ! effects in the internals (e.g. would need to abstract More-Sorensen etc.)
  subroutine mult_inner_dense(A, out_ATA)
    ! Compute A^T * A and store in out_A
    class(dense_matrix), intent(in) :: A
    real(wp), intent(out) :: out_ATA(:,:)
    ! Takes an m x n matrix J and forms the
    ! n x n matrix A given by
    ! A = J' * J

    If (A%fortran_order) then
      call PREC(gemm)('T','N',A%n, A%n, A%m, 1.0_wp, A%data, A%m, A%data, A%m, 0.0_wp, out_ATA, A%n)
    else
      ! c format
      call PREC(gemm)('N','T',A%n, A%n, A%m, 1.0_wp, A%data, A%n, A%data, A%n, 0.0_wp, out_ATA, A%n)
    End If

  end subroutine mult_inner_dense

end module MODULE_PREC(ral_nlls_matrix)
