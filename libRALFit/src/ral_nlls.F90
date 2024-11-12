! Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2020 Numerical Algorithms Group (NAG). All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.


! ral_nlls - RALFit a nonlinear least squares solver

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls)

  use MODULE_PREC(ral_nlls_internal)
  use MODULE_PREC(ral_nlls_workspaces), only : lp, np, wp, params_base_type,           &
    nlls_options, nlls_inform, nlls_workspace, eval_f_type, eval_j_type, eval_hf_type, &
    eval_hp_type
  use MODULE_PREC(ral_nlls_fd), only: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, &
    jacobian_handle, jacobian_setup, jacobian_calc, jacobian_free

  implicit none

  public :: wp ! working precion of the library
  public :: lp, np ! low and normal precion kinds
  ! Abstract interfaces for call-backs
  public :: eval_f_type, eval_j_type, eval_hf_type, eval_hp_type
  public :: nlls_solve, nlls_iterate, nlls_finalize, nlls_strerror
  public :: nlls_options, nlls_inform, nlls_workspace
  public :: params_base_type
  public :: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy
  public :: jacobian_handle, jacobian_setup, jacobian_calc, jacobian_free

end module MODULE_PREC(ral_nlls)

