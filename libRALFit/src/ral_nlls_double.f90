! Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_double :: a nonlinear least squares solver

module ral_nlls_double

  use ral_nlls_internal
  use ral_nlls_workspaces, only : params_base_type,           &
    nlls_options, nlls_inform, nlls_workspace
  use ral_nlls_fd, only: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, &
    jacobian_handle, jacobian_setup, jacobian_calc, jacobian_free

  implicit none

  private

  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: f
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hf_type
  end interface

  public :: nlls_solve, nlls_iterate, nlls_finalize, nlls_strerror
  public :: nlls_options, nlls_inform, nlls_workspace
  public :: params_base_type
  public :: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy
  public :: jacobian_handle, jacobian_setup, jacobian_calc, jacobian_free

end module ral_nlls_double

