module ral_nlls_ciface

  use iso_c_binding
  use ral_nlls_double, only:                       &
       f_nlls_options      => nlls_options,        &
       f_nlls_inform       => nlls_inform,         &
       f_nlls_workspace    => nlls_workspace,      &
       f_nlls_solve        => nlls_solve,          &
       f_nlls_iterate      => nlls_iterate,        &
       f_params_base_type  => params_base_type
  implicit none

  integer, parameter :: wp = C_DOUBLE

  type, bind(C) :: nlls_options
     integer(C_INT) :: f_arrays ! true (!=0) or false (==0)

     integer(C_INT) :: out
     integer(C_INT) :: print_level
     logical(c_bool) :: print_options
     integer(C_INT) :: print_header
     integer(C_INT) :: maxit
     integer(C_INT) :: model
     integer(C_INT) :: type_of_method
     integer(C_INT) :: nlls_method
     logical(c_bool) :: allow_fallback_method
     integer(C_INT) :: lls_solver
     real(wp) :: stop_g_absolute
     real(wp) :: stop_g_relative
     real(wp) :: stop_f_absolute
     real(wp) :: stop_f_relative
     real(wp) :: stop_s
     integer(c_int) :: relative_tr_radius
     real(wp) :: initial_radius_scale
     real(wp) :: initial_radius
     real(wp) :: base_regularization
     integer(c_int) :: regularization
     real(wp) :: regularization_term
     real(wp) :: regularization_power
     real(wp) :: maximum_radius
     real(wp) :: eta_successful
     real(wp) :: eta_success_but_reduce
     real(wp) :: eta_very_successful
     real(wp) :: eta_too_successful
     real(wp) :: radius_increase
     real(wp) :: radius_reduce
     real(wp) :: radius_reduce_max
     integer(c_int) :: tr_update_strategy
     real(wp) :: hybrid_switch
     logical(c_bool) :: exact_second_derivatives
     logical(c_bool) :: subproblem_eig_fact
     logical(c_bool) :: use_ews_subproblem
     logical(c_bool) :: force_min_eig_symm
     integer(C_INT) :: scale
     real(wp) :: scale_max
     real(wp) :: scale_min
     logical(c_bool) :: scale_trim_min
     logical(c_bool) :: scale_trim_max
     logical(c_bool) :: scale_require_increase
     logical(c_bool) :: setup_workspaces
     logical(c_bool) :: remove_workspaces
     integer(c_int) :: more_sorensen_maxits
     real(wp) :: more_sorensen_shift
     real(wp) :: more_sorensen_tiny
     real(wp) :: more_sorensen_tol
     real(wp) :: hybrid_tol
     integer(c_int) :: hybrid_switch_its
     real(wp) :: reg_order
     integer(c_int) :: inner_method
     logical(c_bool) :: output_progress_vectors
     logical(c_bool) :: update_lower_order
     logical(c_bool) :: Fortran_Jacobian

     integer(c_int) :: box_nFref_max
     real(wp) :: box_gamma
     real(wp) :: box_decmin
     real(wp) :: box_bigbnd
     real(wp) :: box_wolfe_descent
     real(wp) :: box_wolfe_curvature
     real(wp) :: box_kanzow_power
     real(wp) :: box_kanzow_descent
     real(wp) :: box_quad_model_descent
     Logical(c_bool):: box_tr_test_step
     Logical(c_bool):: box_wolfe_test_step
     real(wp) :: box_tau_descent
     integer(c_int):: box_max_ntrfail
     integer(c_int):: box_quad_match
     real(wp) :: box_alpha_scale
     real(wp) :: box_Delta_scale
     real(wp) :: box_tau_min
     integer(c_int):: box_ls_step_maxit
     integer(c_int):: box_linesearch_type
  end type nlls_options

  type, bind(C) :: nlls_inform
     integer(C_INT) :: status
     character( kind = c_char), dimension(81) :: error_message
     integer(C_INT) :: alloc_status
     character( kind = c_char), dimension(81) :: bad_alloc
     integer(C_INT) :: iter
     integer(C_INT) :: inner_iter
     LOGICAL(c_bool) :: inner_iter_success
     integer(C_INT) :: f_eval
     integer(C_INT) :: g_eval
     integer(C_INT) :: h_eval
     integer(C_INT) :: convergence_normf
     integer(C_INT) :: convergence_normg
     integer(C_INT) :: convergence_norms
     real(wp) :: resinf
     real(wp) :: gradinf
     real(wp) :: obj
     real(wp) :: norm_g
     real(wp) :: scaled_g
     integer(C_INT) :: external_return
     character( kind = c_char), dimension(81) :: external_name
     real(wp) :: step
     Integer(c_int) :: ls_step_iter
     Integer(c_int) :: f_eval_ls
     Integer(c_int) :: g_eval_ls
     Integer(c_int) :: pg_step_iter
     Integer(c_int) :: f_eval_pg
     Integer(c_int) :: g_eval_pg
  end type nlls_inform

  abstract interface
     integer(c_int) function c_eval_r_type(n, m, params, x, r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: r
     end function c_eval_r_type
  end interface

  abstract interface
     integer(c_int) function c_eval_j_type(n, m, params, x, j) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: j
     end function c_eval_j_type
  end interface

  abstract interface
     integer(c_int) function c_eval_hf_type(n, m, params, x, f, hf) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(in) :: f
       real(c_double), dimension(*), intent(out) :: hf
     end function c_eval_hf_type
  end interface

  type, extends(f_params_base_type) :: params_wrapper
     procedure(c_eval_r_type), nopass, pointer :: r
     procedure(c_eval_j_type), nopass, pointer :: j
     procedure(c_eval_hf_type), nopass, pointer :: hf
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
  coptions%eta_very_successful = foptions%eta_very_successful
  coptions%eta_too_successful = foptions%eta_too_successful
  coptions%radius_increase = foptions%radius_increase
  coptions%radius_reduce = foptions%radius_reduce
  coptions%radius_reduce_max = foptions%radius_reduce_max
  coptions%tr_update_strategy = foptions%tr_update_strategy
  coptions%hybrid_switch = foptions%hybrid_switch
  coptions%exact_second_derivatives = foptions%exact_second_derivatives
  coptions%subproblem_eig_fact = foptions%subproblem_eig_fact
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

end subroutine ral_nlls_default_options_d

subroutine nlls_solve_d(n, m, cx, r, j, hf,  params, coptions, cinform, cweights) bind(C)
  use ral_nlls_ciface
  implicit none

  integer( C_INT ) , INTENT( IN ), value :: n, m
  real( c_double ), dimension(*) :: cx
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
  real( c_double ), dimension(:), pointer :: fweights
!  real( wp ), dimension(*), optional :: cweights

  logical :: f_arrays

  ! copy data in and associate pointers correctly
  call copy_options_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  fparams%params = params

  if (C_ASSOCIATED(cweights)) then
     call c_f_pointer(cweights, fweights, shape = (/ m /) )
     call f_nlls_solve( n, m, cx, &
       c_eval_r, c_eval_j,   &
       c_eval_hf, fparams,   &
       foptions,finform,fweights)
  else
     call f_nlls_solve( n, m, cx, &
       c_eval_r, c_eval_j,   &
       c_eval_hf, fparams,   &
       foptions,finform)
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

   if(c_associated(cw)) return ! Nothing to do

   call c_f_pointer(cw, fw)
   deallocate(fw)
   cw = C_NULL_PTR
end subroutine ral_nlls_free_workspace_d

subroutine ral_nlls_iterate_d(n, m, cx, cw, r, j, hf, params, coptions, &
      cinform, cweights) bind(C)
  use ral_nlls_ciface
  implicit none

  integer( C_INT) , INTENT( IN ), value :: n, m
  real( c_double ), dimension(*) :: cx
  type( C_PTR), value :: cw
  type( C_FUNPTR ), value :: r
  type( C_FUNPTR ), value :: j
  type( C_FUNPTR ), value :: hf
  type( C_PTR ), value :: params
  TYPE( nlls_options ) :: coptions
  TYPE( nlls_inform )  :: cinform
  TYPE( C_PTR ), value :: cweights

  type( params_wrapper ) :: fparams
  TYPE( f_nlls_inform) :: finform
  TYPE( f_nlls_workspace ), pointer :: fw
  real( c_double ), dimension(:), pointer :: fweights
  TYPE( f_nlls_options) :: foptions

  logical :: f_arrays

  ! copy data in and associate pointers correctly
  call copy_options_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  call c_f_pointer(cw, fw)
  fparams%params = params

  if (C_ASSOCIATED(cweights)) then
     call c_f_pointer(cweights, fweights, shape = (/ m /) )
     call f_nlls_iterate( n, m, cx, fw, &
          c_eval_r, c_eval_j,   &
          c_eval_hf, fparams,   &
          finform, foptions, fweights)
  else
     call f_nlls_iterate( n, m, cx, fw, &
          c_eval_r, c_eval_j,   &
          c_eval_hf, fparams,   &
          finform, foptions)
  end if

  ! Copy data out
  call copy_info_out(finform, cinform)

end subroutine ral_nlls_iterate_d
