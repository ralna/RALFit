! Copyright (c) 2015, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
module ral_nlls_ciface

  use ral_nlls_types
  use iso_c_binding
  use ral_nlls_double, only:                       &
       f_nlls_options      => nlls_options,        &
       f_nlls_inform       => nlls_inform,         &
       f_nlls_workspace    => nlls_workspace,      &
       f_nlls_solve        => nlls_solve,          &
       f_nlls_iterate      => nlls_iterate,        &
       f_params_base_type  => params_base_type
  implicit none

  type, bind(C) :: nlls_options
     integer(ral_c_int) :: f_arrays ! true (!=0) or false (==0)

     integer(ral_c_int) :: out
     integer(ral_c_int) :: print_level
     logical(c_bool) :: print_options
     integer(ral_c_int) :: print_header
     integer(ral_c_int) :: maxit
     integer(ral_c_int) :: model
     integer(ral_c_int) :: type_of_method
     integer(ral_c_int) :: nlls_method
     logical(c_bool) :: allow_fallback_method
     integer(ral_c_int) :: lls_solver
     real(ral_c_real) :: stop_g_absolute
     real(ral_c_real) :: stop_g_relative
     real(ral_c_real) :: stop_f_absolute
     real(ral_c_real) :: stop_f_relative
     real(ral_c_real) :: stop_s
     integer(ral_c_int) :: relative_tr_radius
     real(ral_c_real) :: initial_radius_scale
     real(ral_c_real) :: initial_radius
     real(ral_c_real) :: base_regularization
     integer(ral_c_int) :: regularization
     real(ral_c_real) :: regularization_term
     real(ral_c_real) :: regularization_power
     real(ral_c_real) :: maximum_radius
     real(ral_c_real) :: eta_successful
     real(ral_c_real) :: eta_success_but_reduce
     real(ral_c_real) :: eta_very_successful
     real(ral_c_real) :: eta_too_successful
     real(ral_c_real) :: radius_increase
     real(ral_c_real) :: radius_reduce
     real(ral_c_real) :: radius_reduce_max
     integer(ral_c_int) :: tr_update_strategy
     real(ral_c_real) :: hybrid_switch
     logical(c_bool) :: exact_second_derivatives
     logical(c_bool) :: subproblem_eig_fact
     logical(c_bool) :: use_ews_subproblem
     logical(c_bool) :: force_min_eig_symm
     integer(ral_c_int) :: scale
     real(ral_c_real) :: scale_max
     real(ral_c_real) :: scale_min
     logical(c_bool) :: scale_trim_min
     logical(c_bool) :: scale_trim_max
     logical(c_bool) :: scale_require_increase
     logical(c_bool) :: setup_workspaces
     logical(c_bool) :: remove_workspaces
     integer(ral_c_int) :: more_sorensen_maxits
     real(ral_c_real) :: more_sorensen_shift
     real(ral_c_real) :: more_sorensen_tiny
     real(ral_c_real) :: more_sorensen_tol
     real(ral_c_real) :: hybrid_tol
     integer(ral_c_int) :: hybrid_switch_its
     real(ral_c_real) :: reg_order
     integer(ral_c_int) :: inner_method
     logical(c_bool) :: output_progress_vectors
     logical(c_bool) :: update_lower_order
     logical(c_bool) :: Fortran_Jacobian

     integer(ral_c_int) :: box_nFref_max
     real(ral_c_real) :: box_gamma
     real(ral_c_real) :: box_decmin
     real(ral_c_real) :: box_bigbnd
     real(ral_c_real) :: box_wolfe_descent
     real(ral_c_real) :: box_wolfe_curvature
     real(ral_c_real) :: box_kanzow_power
     real(ral_c_real) :: box_kanzow_descent
     real(ral_c_real) :: box_quad_model_descent
     Logical(c_bool):: box_tr_test_step
     Logical(c_bool):: box_wolfe_test_step
     real(ral_c_real) :: box_tau_descent
     integer(ral_c_int):: box_max_ntrfail
     integer(ral_c_int):: box_quad_match
     real(ral_c_real) :: box_alpha_scale
     real(ral_c_real) :: box_Delta_scale
     real(ral_c_real) :: box_tau_min
     integer(ral_c_int):: box_ls_step_maxit
     integer(ral_c_int):: box_linesearch_type
     real(ral_c_real) :: fd_step
     Integer(ral_c_int) :: check_derivatives
     Real(ral_c_real) :: derivative_test_tol
  end type nlls_options

  type, bind(C) :: nlls_inform
     integer(ral_c_int) :: status
     character( kind = c_char), dimension(81) :: error_message
     integer(ral_c_int) :: alloc_status
     character( kind = c_char), dimension(81) :: bad_alloc
     integer(ral_c_int) :: iter
     integer(ral_c_int) :: inner_iter
     LOGICAL(c_bool) :: inner_iter_success
     integer(ral_c_int) :: f_eval
     integer(ral_c_int) :: g_eval
     integer(ral_c_int) :: h_eval
     integer(ral_c_int) :: hp_eval
     integer(ral_c_int) :: convergence_normf
     integer(ral_c_int) :: convergence_normg
     integer(ral_c_int) :: convergence_norms
     real(ral_c_real) :: resinf
     real(ral_c_real) :: gradinf
     real(ral_c_real) :: obj
     real(ral_c_real) :: norm_g
     real(ral_c_real) :: scaled_g
     integer(ral_c_int) :: external_return
     character( kind = c_char), dimension(81) :: external_name
     real(ral_c_real) :: step
     Integer(ral_c_int) :: ls_step_iter
     Integer(ral_c_int) :: f_eval_ls
     Integer(ral_c_int) :: g_eval_ls
     Integer(ral_c_int) :: pg_step_iter
     Integer(ral_c_int) :: f_eval_pg
     Integer(ral_c_int) :: g_eval_pg
     Integer(ral_c_int) :: fd_f_eval
  end type nlls_inform

  abstract interface
     integer(ral_c_int) function c_eval_r_type(n, m, params, x, r) bind(c)
       use, intrinsic :: iso_c_binding
       use :: ral_nlls_types
       implicit none
       integer(ral_c_int), value :: n, m
       type(C_PTR), value :: params
       real(ral_c_real), dimension(*), intent(in) :: x
       real(ral_c_real), dimension(*), intent(out) :: r
     end function c_eval_r_type
  end interface

  abstract interface
     integer(ral_c_int) function c_eval_j_type(n, m, params, x, j) bind(c)
       use, intrinsic :: iso_c_binding
       use :: ral_nlls_types
       implicit none
       integer(ral_c_int), value :: n,m
       type(C_PTR), value :: params
       real(ral_c_real), dimension(*), intent(in) :: x
       real(ral_c_real), dimension(*), intent(out) :: j
     end function c_eval_j_type
  end interface

  abstract interface
     integer(ral_c_int) function c_eval_hf_type(n, m, params, x, f, hf) bind(c)
       use, intrinsic :: iso_c_binding
       use :: ral_nlls_types
       implicit none
       integer(ral_c_int), value :: n,m
       type(C_PTR), value :: params
       real(ral_c_real), dimension(*), intent(in) :: x
       real(ral_c_real), dimension(*), intent(in) :: f
       real(ral_c_real), dimension(*), intent(out) :: hf
     end function c_eval_hf_type
  end interface

  abstract interface
     integer(ral_c_int) function c_eval_hp_type(n, m, x, y, hp, params) bind(c)
       use, intrinsic :: iso_c_binding
       use :: ral_nlls_types
       implicit none
       integer(ral_c_int) :: n,m
       real(ral_c_real), dimension(*), intent(in)  :: x
       real(ral_c_real), dimension(*), intent(in)  :: y
       real(ral_c_real), dimension(*), intent(out) :: hp
       type(C_PTR), value :: params
     end function c_eval_hp_type
  end interface

  type, extends(f_params_base_type) :: params_wrapper
     procedure(c_eval_r_type), nopass, pointer :: r
     procedure(c_eval_j_type), nopass, pointer :: j
     procedure(c_eval_hf_type), nopass, pointer :: hf
     procedure(c_eval_hp_type), nopass, pointer :: hp
     type(C_PTR) :: params
  end type params_wrapper

contains


  subroutine copy_options_in(coptions, foptions, f_arrays)

    type( nlls_options ), intent(in) :: coptions
    type( f_nlls_options ), intent(out) :: foptions
    logical, intent(out) :: f_arrays

    f_arrays = (coptions%f_arrays .ne. 0)
    foptions%out = coptions%out
    foptions%print_level = coptions%print_level
    foptions%print_options = coptions%print_options
    foptions%print_header = coptions%print_header
    foptions%maxit = coptions%maxit
    foptions%model = coptions%model
    foptions%type_of_method = coptions%type_of_method
    foptions%nlls_method = coptions%nlls_method
    foptions%allow_fallback_method = coptions%allow_fallback_method
    foptions%lls_solver = coptions%lls_solver
    foptions%stop_g_absolute = coptions%stop_g_absolute
    foptions%stop_g_relative = coptions%stop_g_relative
    foptions%stop_f_absolute = coptions%stop_f_absolute
    foptions%stop_f_relative = coptions%stop_f_relative
    foptions%stop_s = coptions%stop_s
    foptions%relative_tr_radius = coptions%relative_tr_radius
    foptions%initial_radius_scale = coptions%initial_radius_scale
    foptions%initial_radius = coptions%initial_radius
    foptions%base_regularization = coptions%base_regularization
    foptions%regularization = coptions%regularization
    foptions%regularization_term = coptions%regularization_term
    foptions%regularization_power = coptions%regularization_power
    foptions%maximum_radius = coptions%maximum_radius
    foptions%eta_successful = coptions%eta_successful
    foptions%eta_success_but_reduce = coptions%eta_success_but_reduce
    foptions%eta_very_successful = coptions%eta_very_successful
    foptions%eta_too_successful = coptions%eta_too_successful
    foptions%radius_increase = coptions%radius_increase
    foptions%radius_reduce = coptions%radius_reduce
    foptions%radius_reduce_max = coptions%radius_reduce_max
    foptions%tr_update_strategy = coptions%tr_update_strategy
    foptions%hybrid_switch = coptions%hybrid_switch
    foptions%exact_second_derivatives = coptions%exact_second_derivatives
    foptions%subproblem_eig_fact = coptions%subproblem_eig_fact
    foptions%use_ews_subproblem = coptions%use_ews_subproblem
    foptions%force_min_eig_symm = coptions%force_min_eig_symm
    foptions%scale = coptions%scale
    foptions%scale_max = coptions%scale_max
    foptions%scale_min = coptions%scale_min
    foptions%scale_trim_max = coptions%scale_trim_max
    foptions%scale_trim_min = coptions%scale_trim_min
    foptions%scale_require_increase = coptions%scale_require_increase
    foptions%setup_workspaces = coptions%setup_workspaces
    foptions%remove_workspaces = coptions%remove_workspaces
    foptions%more_sorensen_maxits = coptions%more_sorensen_maxits
    foptions%more_sorensen_shift = coptions%more_sorensen_shift
    foptions%more_sorensen_tiny = coptions%more_sorensen_tiny
    foptions%more_sorensen_tol = coptions%more_sorensen_tol
    foptions%hybrid_tol = coptions%hybrid_tol
    foptions%hybrid_switch_its = coptions%hybrid_switch_its
    foptions%reg_order = coptions%reg_order
    foptions%inner_method = coptions%inner_method
    foptions%output_progress_vectors = coptions%output_progress_vectors
    foptions%update_lower_order = coptions%update_lower_order
    foptions%Fortran_Jacobian = coptions%Fortran_Jacobian
    foptions%box_nFref_max = coptions%box_nFref_max
    foptions%box_gamma = coptions%box_gamma
    foptions%box_decmin = coptions%box_decmin
    foptions%box_bigbnd = coptions%box_bigbnd
    foptions%box_wolfe_descent = coptions%box_wolfe_descent
    foptions%box_wolfe_curvature = coptions%box_wolfe_curvature
    foptions%box_kanzow_power = coptions%box_kanzow_power
    foptions%box_kanzow_descent = coptions%box_kanzow_descent
    foptions%box_quad_model_descent = coptions%box_quad_model_descent
    foptions%box_tr_test_step = coptions%box_tr_test_step
    foptions%box_wolfe_test_step = coptions%box_wolfe_test_step
    foptions%box_tau_descent = coptions%box_tau_descent
    foptions%box_max_ntrfail = coptions%box_max_ntrfail
    foptions%box_quad_match = coptions%box_quad_match
    foptions%box_alpha_scale = coptions%box_alpha_scale
    foptions%box_Delta_scale = coptions%box_Delta_scale
    foptions%box_tau_min = coptions%box_tau_min
    foptions%box_ls_step_maxit = coptions%box_ls_step_maxit
    foptions%box_linesearch_type = coptions%box_linesearch_type
    foptions%fd_step = coptions%fd_step
    foptions%check_derivatives = coptions%check_derivatives
    foptions%derivative_test_tol = coptions%derivative_test_tol

  end subroutine copy_options_in

  subroutine copy_info_out(finfo,cinfo)

    type(f_nlls_inform), intent(in) :: finfo
    type(nlls_inform) , intent(out) :: cinfo

    integer :: i

    cinfo%status = finfo%status
    do i = 1,len(finfo%error_message)
       cinfo%error_message(i) = finfo%error_message(i:i)
    end do
    cinfo%error_message(len(finfo%error_message) + 1) = C_NULL_CHAR
    cinfo%alloc_status = finfo%alloc_status
    do i = 1,len(finfo%bad_alloc)
       cinfo%bad_alloc(i) = finfo%bad_alloc(i:i)
    end do
    cinfo%bad_alloc(len(finfo%bad_alloc) + 1) = C_NULL_CHAR
    cinfo%iter = finfo%iter
    cinfo%step = finfo%step
    cinfo%inner_iter = finfo%iter
    cinfo%f_eval = finfo%f_eval
    cinfo%g_eval = finfo%g_eval
    cinfo%h_eval = finfo%h_eval
    cinfo%hp_eval = finfo%hp_eval
    cinfo%convergence_normf = finfo%convergence_normf
    cinfo%convergence_normg = finfo%convergence_normg
    cinfo%convergence_norms = finfo%convergence_norms
    if(allocated(finfo%resvec)) &
       cinfo%resinf = maxval(abs(finfo%resvec(:)))
    if(allocated(finfo%gradvec)) &
       cinfo%gradinf = maxval(abs(finfo%gradvec(:)))
    cinfo%obj = finfo%obj
    cinfo%norm_g = finfo%norm_g
    cinfo%scaled_g = finfo%scaled_g
    cinfo%external_return = finfo%external_return
    do i = 1,len(finfo%external_name)
       cinfo%external_name(i) = finfo%external_name(i:i)
    end do
    cinfo%external_name(len(finfo%external_name) + 1) = C_NULL_CHAR
    cinfo%step = finfo%step
    cinfo%ls_step_iter = finfo%ls_step_iter
    cinfo%f_eval_ls = finfo%f_eval_ls
    cinfo%g_eval_ls = finfo%g_eval_ls
    cinfo%pg_step_iter = finfo%pg_step_iter
    cinfo%f_eval_pg = finfo%f_eval_pg
    cinfo%g_eval_pg = finfo%g_eval_pg
    cinfo%fd_f_eval = finfo%fd_f_eval

  end subroutine copy_info_out

  subroutine c_eval_r(evalrstatus, n, m, x, f, fparams)
    integer, intent(in) :: n, m
    integer, intent(out) :: evalrstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(out) :: f
    class(f_params_base_type), intent(inout) :: fparams

    select type(fparams)
    type is(params_wrapper)
       evalrstatus =  fparams%r(n,m,fparams%params,x(1:n),f(1:m))
    end select

  end subroutine c_eval_r

  subroutine c_eval_j(evaljstatus, n, m, x, j, fparams)
    integer, intent(in) :: n, m
    integer, intent(out) :: evaljstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(out) :: j
    class(f_params_base_type), intent(inout) :: fparams

    select type(fparams)
    type is(params_wrapper)
       evaljstatus = fparams%j(n,m,fparams%params,x(1:n),j(1:n*m))
    end select

  end subroutine c_eval_j

  subroutine c_eval_hf(evalhstatus, n, m, x, f, hf, fparams)
    integer, intent(in) :: n, m
    integer, intent(out) :: evalhstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(in) :: f
    double precision, dimension(*), intent(out) :: hf
    class(f_params_base_type), intent(inout) :: fparams

    select type(fparams)
    type is(params_wrapper)
       evalhstatus = fparams%hf(n,m,fparams%params,x(1:n),f(1:m),hf(1:n**2))
    end select

  end subroutine c_eval_hf

  subroutine c_eval_hp(evalhpstatus, n, m, x, y, hp, fparams)
    integer, intent(in) :: n, m
    integer, intent(out) :: evalhpstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(in) :: y
    double precision, dimension(*), intent(out) :: hp
    class(f_params_base_type), intent(inout) :: fparams

    select type(fparams)
    type is(params_wrapper)
       evalhpstatus = fparams%hp(n,m,x(1:n),y(1:n),hp(1:n*m),fparams%params)
    end select

  end subroutine c_eval_hp

  integer(c_int) function c_eval_j_dummy(n, m, params, x, j) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: j

       continue
       c_eval_j_dummy = -45544554 ! Magic number to request FD
  end function c_eval_j_dummy

  integer(c_int) function c_eval_hf_dummy(n, m, params, x, f, hf) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(in) :: f
       real(c_double), dimension(*), intent(out) :: hf
       continue
       c_eval_hf_dummy =  -1023
  end function c_eval_hf_dummy

end module ral_nlls_ciface

subroutine ral_nlls_default_options_d(coptions) bind(C)
  use ral_nlls_ciface
  implicit none

  type(nlls_options), intent(out) :: coptions
  type(f_nlls_options) :: foptions


  coptions%f_arrays = 0 ! (false) default to C style arrays
  coptions%out = foptions%out
  coptions%print_level = foptions%print_level
  coptions%print_options = foptions%print_options
  coptions%print_header = foptions%print_header
  coptions%maxit = foptions%maxit
  coptions%model = foptions%model
  coptions%type_of_method = foptions%type_of_method
  coptions%nlls_method = foptions%nlls_method
  coptions%allow_fallback_method = foptions%allow_fallback_method
  coptions%lls_solver = foptions%lls_solver
  coptions%stop_g_absolute = foptions%stop_g_absolute
  coptions%stop_g_relative = foptions%stop_g_relative
  coptions%stop_f_absolute = foptions%stop_f_absolute
  coptions%stop_f_relative = foptions%stop_f_relative
  coptions%stop_s = foptions%stop_s
  coptions%relative_tr_radius = foptions%relative_tr_radius
  coptions%initial_radius_scale = foptions%initial_radius_scale
  coptions%initial_radius = foptions%initial_radius
  coptions%base_regularization = foptions%base_regularization
  coptions%regularization = foptions%regularization
  coptions%regularization_term = foptions%regularization_term
  coptions%regularization_power = foptions%regularization_power
  coptions%maximum_radius = foptions%maximum_radius
  coptions%eta_successful = foptions%eta_successful
  coptions%eta_success_but_reduce = foptions%eta_success_but_reduce
  coptions%eta_very_successful = foptions%eta_very_successful
  coptions%eta_too_successful = foptions%eta_too_successful
  coptions%radius_increase = foptions%radius_increase
  coptions%radius_reduce = foptions%radius_reduce
  coptions%radius_reduce_max = foptions%radius_reduce_max
  coptions%tr_update_strategy = foptions%tr_update_strategy
  coptions%hybrid_switch = foptions%hybrid_switch
  coptions%exact_second_derivatives = foptions%exact_second_derivatives
  coptions%subproblem_eig_fact = foptions%subproblem_eig_fact
  coptions%use_ews_subproblem = foptions%use_ews_subproblem
  coptions%force_min_eig_symm = foptions%force_min_eig_symm
  coptions%scale = foptions%scale
  coptions%scale_max = foptions%scale_max
  coptions%scale_min = foptions%scale_min
  coptions%scale_trim_max = foptions%scale_trim_max
  coptions%scale_trim_min = foptions%scale_trim_min
  coptions%scale_require_increase = foptions%scale_require_increase
  coptions%setup_workspaces = foptions%setup_workspaces
  coptions%remove_workspaces = foptions%remove_workspaces
  coptions%more_sorensen_maxits = foptions%more_sorensen_maxits
  coptions%more_sorensen_shift = foptions%more_sorensen_shift
  coptions%more_sorensen_tiny = foptions%more_sorensen_tiny
  coptions%more_sorensen_tol = foptions%more_sorensen_tol
  coptions%hybrid_tol = foptions%hybrid_tol
  coptions%hybrid_switch_its = foptions%hybrid_switch_its
  coptions%reg_order = foptions%reg_order
  coptions%inner_method = foptions%inner_method
  coptions%output_progress_vectors = foptions%output_progress_vectors
  coptions%update_lower_order = foptions%update_lower_order
  coptions%Fortran_Jacobian = foptions%Fortran_Jacobian

  coptions%box_nFref_max = foptions%box_nFref_max
  coptions%box_gamma = foptions%box_gamma
  coptions%box_decmin = foptions%box_decmin
  coptions%box_bigbnd = foptions%box_bigbnd
  coptions%box_wolfe_descent = foptions%box_wolfe_descent
  coptions%box_wolfe_curvature = foptions%box_wolfe_curvature
  coptions%box_kanzow_power = foptions%box_kanzow_power
  coptions%box_kanzow_descent = foptions%box_kanzow_descent
  coptions%box_quad_model_descent = foptions%box_quad_model_descent
  coptions%box_tr_test_step = foptions%box_tr_test_step
  coptions%box_wolfe_test_step = foptions%box_wolfe_test_step
  coptions%box_tau_descent = foptions%box_tau_descent
  coptions%box_max_ntrfail = foptions%box_max_ntrfail
  coptions%box_quad_match = foptions%box_quad_match
  coptions%box_alpha_scale = foptions%box_alpha_scale
  coptions%box_Delta_scale = foptions%box_Delta_scale
  coptions%box_tau_min = foptions%box_tau_min
  coptions%box_ls_step_maxit = foptions%box_ls_step_maxit
  coptions%box_linesearch_type = foptions%box_linesearch_type

  coptions%fd_step = foptions%fd_step
  coptions%check_derivatives = foptions%check_derivatives
  coptions%derivative_test_tol = foptions%derivative_test_tol

end subroutine ral_nlls_default_options_d

subroutine nlls_solve_d(n, m, cx, r, j, hf,  params, coptions, cinform, &
     cweights, hp, clower_bounds, cupper_bounds) bind(C)
  use ral_nlls_ciface
  use ral_nlls_double, only: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy
  implicit none

  integer( ral_c_int ) , INTENT( IN ), value :: n, m
  real( ral_c_real ), dimension(*) :: cx
  type( C_FUNPTR ), value :: r
  type( C_FUNPTR ), value :: j
  type( C_FUNPTR ), value :: hf
  type( C_PTR ), value :: params
  TYPE( nlls_inform )  :: cinform
  TYPE( nlls_options ) :: coptions
  type( params_wrapper ) :: fparams
  TYPE( f_nlls_options ) :: foptions
  TYPE( f_nlls_inform ) :: finform
  TYPE( C_PTR ), value :: cweights
  TYPE( C_PTR ), value :: clower_bounds
  TYPE( C_PTR ), value :: cupper_bounds
  real( ral_c_real ), dimension(:), pointer :: fweights
  type( C_FUNPTR ), value :: hp
  real( ral_c_real ), dimension(:), pointer :: flower_bounds
  real( ral_c_real ), dimension(:), pointer :: fupper_bounds

!  real( wp ), dimension(*), optional :: cweights

  logical :: f_arrays
  logical :: lb_sent_in, ub_sent_in
  logical :: hp_sent_in = .false.

  ! copy data in and associate pointers correctly
  call copy_options_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)

  if (c_associated(j)) then
    call c_f_procpointer(j, fparams%j)
  else
    fparams%j => c_eval_j_dummy
  end if

  if (c_associated(hf)) then
    call c_f_procpointer(hf, fparams%hf)
  else
    fparams%hf => c_eval_hf_dummy
  endif

  fparams%params = params
  if (C_ASSOCIATED(hp)) hp_sent_in = .true.

  ! the following steps for passing optional arguments
  ! requires a compiler compatible with Fortran 2008+TS29113
  nullify(fweights)
  if (C_ASSOCIATED(cweights)) then
     call c_f_pointer(cweights, fweights, shape = (/ m /) )
  end if
  ! In the Fortran standard, passing null pointers as an optional argument is valid
  ! and is considered not present.
  ! However, there seems to be a bug in aocc5.0 where that behaviour is inconsistent.
  ! A lot of the below if statements can be replaced once this is fixed.
  nullify(flower_bounds)
  lb_sent_in = .false.
  if (C_ASSOCIATED(clower_bounds)) then
     lb_sent_in = .true.
     call c_f_pointer(clower_bounds, flower_bounds, shape = (/ n /) )
  end if
  nullify(fupper_bounds)
  ub_sent_in = .false.
  if (C_ASSOCIATED(cupper_bounds)) then
     ub_sent_in = .true.
     call c_f_pointer(cupper_bounds, fupper_bounds, shape = (/ n /) )
  end if


  if (hp_sent_in) then
      if (lb_sent_in .and. ub_sent_in) then
         call f_nlls_solve( n, m, cx, &
               c_eval_r, c_eval_j,   &
               c_eval_hf, fparams,   &
               foptions,finform, &
               weights=fweights, &
               eval_hp=c_eval_hp, &
               lower_bounds=flower_bounds, &
               upper_bounds=fupper_bounds)
      elseif (lb_sent_in) then
         call f_nlls_solve( n, m, cx, &
               c_eval_r, c_eval_j,   &
               c_eval_hf, fparams,   &
               foptions,finform, &
               weights=fweights, &
               eval_hp=c_eval_hp, &
               lower_bounds=flower_bounds)
      elseif (ub_sent_in) then
         call f_nlls_solve( n, m, cx, &
               c_eval_r, c_eval_j,   &
               c_eval_hf, fparams,   &
               foptions,finform, &
               weights=fweights, &
               eval_hp=c_eval_hp, &
               upper_bounds=fupper_bounds)
      else
         call f_nlls_solve( n, m, cx, &
               c_eval_r, c_eval_j,   &
               c_eval_hf, fparams,   &
               foptions,finform, &
               weights=fweights, &
               eval_hp=c_eval_hp)
      endif
  else
     if (lb_sent_in .and. ub_sent_in) then
         call f_nlls_solve( n, m, cx, &
            c_eval_r, c_eval_j,   &
            c_eval_hf, fparams,   &
            foptions,finform, &
            weights=fweights, &
            lower_bounds=flower_bounds, &
            upper_bounds=fupper_bounds)
      elseif (lb_sent_in) then
         call f_nlls_solve( n, m, cx, &
            c_eval_r, c_eval_j,   &
            c_eval_hf, fparams,   &
            foptions,finform, &
            weights=fweights, &
            lower_bounds=flower_bounds)
      elseif (ub_sent_in) then
         call f_nlls_solve( n, m, cx, &
            c_eval_r, c_eval_j,   &
            c_eval_hf, fparams,   &
            foptions,finform, &
            weights=fweights, &
            upper_bounds=fupper_bounds)
      else
         call f_nlls_solve( n, m, cx, &
            c_eval_r, c_eval_j,   &
            c_eval_hf, fparams,   &
            foptions,finform, &
            weights=fweights)
      end if
  end if

  ! Copy data out
   call copy_info_out(finform, cinform)

end subroutine nlls_solve_d

subroutine ral_nlls_init_workspace_d(cw, ciw) bind(C)
   use ral_nlls_ciface
   implicit none

   type(c_ptr) :: cw, ciw

   type(f_nlls_workspace), pointer :: fw, fiw

   allocate(fw)
   allocate(fiw)
!  Link workspace with inner_workspace
   fw%iw_ptr => fiw
!  Self reference inner workspace so not to break recursiveness
   fiw%iw_ptr => fiw

   cw = c_loc(fw)
   ciw = c_loc(fiw)
end subroutine ral_nlls_init_workspace_d

subroutine ral_nlls_free_workspace_d(cw) bind(C)
   use ral_nlls_ciface
   implicit none

   type(c_ptr) :: cw

   type(f_nlls_workspace), pointer :: fw

   if(.not. c_associated(cw)) return ! Nothing to do

   call c_f_pointer(cw, fw)
   deallocate(fw)
   cw = C_NULL_PTR
end subroutine ral_nlls_free_workspace_d

subroutine ral_nlls_iterate_d(n, m, cx, cw, r, j, hf, params, coptions, &
      cinform, cweights, hp, clower_bounds, cupper_bounds) bind(C)
  use ral_nlls_ciface
  implicit none

  integer( ral_c_int) , INTENT( IN ), value :: n, m
  real( ral_c_real ), dimension(*) :: cx
  type( C_PTR), value :: cw
  type( C_FUNPTR ), value :: r
  type( C_FUNPTR ), value :: j
  type( C_FUNPTR ), value :: hf
  type( C_PTR ), value :: params
  TYPE( nlls_options ) :: coptions
  TYPE( nlls_inform )  :: cinform
  TYPE( C_PTR ), value :: cweights
  type( C_FUNPTR ), value :: hp
  type( C_PTR ), value :: clower_bounds
  type( C_PTR ), value :: cupper_bounds

  type( params_wrapper ) :: fparams
  TYPE( f_nlls_inform) :: finform
  TYPE( f_nlls_workspace ), pointer :: fw
  real( ral_c_real ), dimension(:), pointer :: fweights
  real( ral_c_real ), dimension(:), pointer :: flower_bounds
  real( ral_c_real ), dimension(:), pointer :: fupper_bounds
  TYPE( f_nlls_options) :: foptions

  logical :: f_arrays

  ! copy data in and associate pointers correctly
  call copy_options_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  nullify(fparams%hp)
  if (C_associated(hp)) call c_f_procpointer(hp, fparams%hp)
  call c_f_pointer(cw, fw)
  fparams%params = params

  nullify(fweights)
  if (C_ASSOCIATED(cweights)) then
     call c_f_pointer(cweights, fweights, shape = (/ m /) )
  end if
  nullify(flower_bounds)
  if (C_ASSOCIATED(clower_bounds)) then
     call c_f_pointer(clower_bounds, flower_bounds, shape = (/ n /) )
  end if
  nullify(fupper_bounds)
  if (C_ASSOCIATED(cupper_bounds)) then
     call c_f_pointer(cupper_bounds, fupper_bounds, shape = (/ n /) )
  end if

  call f_nlls_iterate( n, m, cx, fw, &
       c_eval_r, c_eval_j,           &
       c_eval_hf, fparams,           &
       finform, foptions,            &
       weights=fweights,             &
       eval_hp=c_eval_hp,            &
       lower_bounds=flower_bounds,   &
       upper_bounds=fupper_bounds)

  ! Copy data out
  call copy_info_out(finform, cinform)

end subroutine ral_nlls_iterate_d
