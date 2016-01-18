module ral_nlls_ciface

  use iso_c_binding
  use ral_nlls_double, only:          &
       f_nlls_control_type => nlls_control_type, &
       f_nlls_inform_type  => nlls_inform_type, &
       f_nlls_workspace    => nlls_workspace, &
       f_ral_nlls => ral_nlls, &
       f_ral_nlls_iterate => ral_nlls_iterate, &
       f_params_base_type => params_base_type
  implicit none

  integer, parameter :: wp = C_DOUBLE

  type, bind(C) :: nlls_control_type
     integer(C_INT) :: f_arrays ! true (!=0) or false (==0)

     integer(C_INT) :: error
     integer(C_INT) :: out
     integer(C_INT) :: print_level
     integer(C_INT) :: maxit
     integer(C_INT) :: model
     integer(C_INT) :: nlls_method
     integer(C_INT) :: lls_solver
     real(wp) :: stop_g_absolute
     real(wp) :: stop_g_relative
     real(wp) :: relative_tr_radius
     real(wp) :: initial_radius_scale
     real(wp) :: initial_radius
     real(wp) :: maximum_radius
     real(wp) :: eta_successful
     real(wp) :: eta_very_successful
     real(wp) :: eta_too_successful
     real(wp) :: radius_increase
     real(wp) :: radius_reduce
     real(wp) :: radius_reduce_max
     real(wp) :: hybrid_switch
     logical(c_bool) :: exact_second_derivatives
     logical(c_bool) :: subproblem_eig_fact
     integer(c_int) :: more_sorensen_maxits
     real(wp) :: more_sorensen_shift
     real(wp) :: more_sorensen_tiny
     real(wp) :: more_sorensen_tol
     real(wp) :: hybrid_tol
     integer(c_int) :: hybrid_switch_its
     logical(c_bool) :: output_progress_vectors
  end type nlls_control_type

  type, bind(C) :: nlls_inform_type
     integer(C_INT) :: status
     integer(C_INT) :: alloc_status
     integer(C_INT) :: iter
     integer(C_INT) :: f_eval
     integer(C_INT) :: g_eval
     integer(C_INT) :: h_eval
     integer(C_INT) :: convergence_normf
     real(wp) :: resinf
     real(wp) :: gradinf
     real(wp) :: obj
     real(wp) :: norm_g
  end type nlls_inform_type

  abstract interface
     integer(c_int) function c_eval_r_type(n, m, params, x, r)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: r
     end function c_eval_r_type
  end interface

  abstract interface
     integer(c_int) function c_eval_j_type(n, m, params, x, j)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       type(C_PTR), value :: params
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: j
     end function c_eval_j_type
  end interface

  abstract interface
     integer(c_int) function c_eval_hf_type(n, m, params, x, f, hf)
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


  subroutine copy_control_in(ccontrol, fcontrol, f_arrays);

    type( nlls_control_type ), intent(in) :: ccontrol
    type( f_nlls_control_type ), intent(out) :: fcontrol
    logical, intent(out) :: f_arrays

    f_arrays = (ccontrol%f_arrays .ne. 0)
    fcontrol%error = ccontrol%error
    fcontrol%out = ccontrol%out
    fcontrol%print_level = ccontrol%print_level
    fcontrol%maxit = ccontrol%maxit
    fcontrol%model = ccontrol%model
    fcontrol%nlls_method = ccontrol%nlls_method
    fcontrol%lls_solver = ccontrol%lls_solver
    fcontrol%stop_g_absolute = ccontrol%stop_g_absolute
    fcontrol%stop_g_relative = ccontrol%stop_g_relative
    fcontrol%initial_radius_scale = ccontrol%initial_radius_scale
    fcontrol%initial_radius = ccontrol%initial_radius
    fcontrol%maximum_radius = ccontrol%maximum_radius
    fcontrol%eta_successful = ccontrol%eta_successful
    fcontrol%eta_very_successful = ccontrol%eta_very_successful
    fcontrol%eta_too_successful = ccontrol%eta_too_successful
    fcontrol%radius_increase = ccontrol%radius_increase
    fcontrol%radius_reduce = ccontrol%radius_reduce
    fcontrol%radius_reduce_max = ccontrol%radius_reduce_max
    fcontrol%hybrid_switch = ccontrol%hybrid_switch
    fcontrol%exact_second_derivatives = ccontrol%exact_second_derivatives
    fcontrol%subproblem_eig_fact = ccontrol%subproblem_eig_fact
    fcontrol%more_sorensen_maxits = ccontrol%more_sorensen_maxits
    fcontrol%more_sorensen_shift = ccontrol%more_sorensen_shift
    fcontrol%more_sorensen_tiny = ccontrol%more_sorensen_tiny
    fcontrol%more_sorensen_tol = ccontrol%more_sorensen_tol
    fcontrol%hybrid_tol = ccontrol%hybrid_tol
    fcontrol%hybrid_switch_its = ccontrol%hybrid_switch_its
    fcontrol%output_progress_vectors = ccontrol%output_progress_vectors

  end subroutine copy_control_in

  subroutine copy_info_out(finfo,cinfo)

    type(f_nlls_inform_type), intent(in) :: finfo
    type(nlls_inform_type) , intent(out) :: cinfo

    cinfo%status = finfo%status
    cinfo%alloc_status = finfo%alloc_status
    cinfo%iter = finfo%iter
    cinfo%f_eval = finfo%f_eval
    cinfo%g_eval = finfo%g_eval
    cinfo%h_eval = finfo%h_eval
    cinfo%convergence_normf = finfo%convergence_normf
    if(allocated(finfo%resvec)) &
       cinfo%resinf = maxval(abs(finfo%resvec(:)))
    if(allocated(finfo%gradvec)) &
       cinfo%gradinf = maxval(abs(finfo%gradvec(:)))
    cinfo%obj = finfo%obj
    cinfo%norm_g = finfo%norm_g

  end subroutine copy_info_out

  subroutine c_eval_r(evalrstatus, n, m, x, f, fparams)
    integer, intent(in) :: n, m
    integer, intent(out) :: evalrstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(out) :: f
    class(f_params_base_type), intent(in) :: fparams

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
    class(f_params_base_type), intent(in) :: fparams

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
    class(f_params_base_type), intent(in) :: fparams

    select type(fparams)
    type is(params_wrapper)
       evalhstatus = fparams%hf(n,m,fparams%params,x(1:n),f(1:m),hf(1:n**2))
    end select

  end subroutine c_eval_hf

end module ral_nlls_ciface

subroutine ral_nlls_default_control_d(ccontrol) bind(C)
  use ral_nlls_ciface
  implicit none

  type(nlls_control_type), intent(out) :: ccontrol
  type(f_nlls_control_type) :: fcontrol


  ccontrol%f_arrays = 0 ! (false) default to C style arrays
  ccontrol%error = fcontrol%error
  ccontrol%out = fcontrol%out
  ccontrol%print_level = fcontrol%print_level
  ccontrol%maxit = fcontrol%maxit
  ccontrol%model = fcontrol%model
  ccontrol%nlls_method = fcontrol%nlls_method
  ccontrol%lls_solver = fcontrol%lls_solver
  ccontrol%stop_g_absolute = fcontrol%stop_g_absolute
  ccontrol%stop_g_relative = fcontrol%stop_g_relative
  ccontrol%relative_tr_radius = fcontrol%relative_tr_radius
  ccontrol%initial_radius_scale = fcontrol%initial_radius_scale
  ccontrol%initial_radius = fcontrol%initial_radius
  ccontrol%maximum_radius = fcontrol%maximum_radius
  ccontrol%eta_successful = fcontrol%eta_successful
  ccontrol%eta_very_successful = fcontrol%eta_very_successful
  ccontrol%eta_too_successful = fcontrol%eta_too_successful
  ccontrol%radius_increase = fcontrol%radius_increase
  ccontrol%radius_reduce = fcontrol%radius_reduce
  ccontrol%radius_reduce_max = fcontrol%radius_reduce_max
  ccontrol%hybrid_switch = fcontrol%hybrid_switch
  ccontrol%exact_second_derivatives = fcontrol%exact_second_derivatives
  ccontrol%subproblem_eig_fact = fcontrol%subproblem_eig_fact
  ccontrol%more_sorensen_maxits = fcontrol%more_sorensen_maxits
  ccontrol%more_sorensen_shift = fcontrol%more_sorensen_shift
  ccontrol%more_sorensen_tiny = fcontrol%more_sorensen_tiny
  ccontrol%more_sorensen_tol = fcontrol%more_sorensen_tol
  ccontrol%hybrid_tol = fcontrol%hybrid_tol
  ccontrol%hybrid_switch_its = fcontrol%hybrid_switch_its
  ccontrol%output_progress_vectors = fcontrol%output_progress_vectors
end subroutine ral_nlls_default_control_d

subroutine ral_nlls_d(n, m, cx, r, j, hf,  params, cinform, coptions) bind(C)
  use ral_nlls_ciface
  implicit none

  integer( C_INT ) , INTENT( IN ), value :: n, m
  real( wp ), dimension(*) :: cx
  type( C_FUNPTR ), value :: r
  type( C_FUNPTR ), value :: j
  type( C_FUNPTR ), value :: hf
  type( C_PTR ), value :: params
  TYPE( nlls_inform_type )  :: cinform
  TYPE( nlls_control_type ) :: coptions

  type( params_wrapper ) :: fparams
  TYPE( f_nlls_inform_type ) :: finform
  TYPE( f_nlls_control_type ) :: foptions

  logical :: f_arrays

  ! copy data in and associate pointers correctly
  call copy_control_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  fparams%params = params

  call f_ral_nlls( n, m, cx, &
       c_eval_r, c_eval_j,   &
       c_eval_hf, fparams,   &
       finform, foptions)

  ! Copy data out
   call copy_info_out(finform, cinform)

end subroutine ral_nlls_d

subroutine ral_nlls_init_workspace_d(cw)
   use ral_nlls_ciface
   implicit none

   type(c_ptr) :: cw

   type(f_nlls_workspace), pointer :: fw

   allocate(fw)
   cw = c_loc(fw)
end subroutine ral_nlls_init_workspace_d

subroutine ral_nlls_free_workspace_d(cw)
   use ral_nlls_ciface
   implicit none

   type(c_ptr) :: cw

   type(f_nlls_workspace), pointer :: fw

   if(c_associated(cw)) return ! Nothing to do

   call c_f_pointer(cw, fw)
   deallocate(fw)
   cw = C_NULL_PTR
end subroutine ral_nlls_free_workspace_d

subroutine ral_nlls_iterate_d(n, m, cx, r, j, hf, params, cinform, coptions, &
      cw) bind(C)
  use ral_nlls_ciface
  implicit none

  integer( C_INT) , INTENT( IN ), value :: n, m
  real( wp ), dimension(*) :: cx
  type( C_FUNPTR ), value :: r
  type( C_FUNPTR ), value :: j
  type( C_FUNPTR ), value :: hf
  type( C_PTR ), value :: params
  TYPE( nlls_inform_type )  :: cinform
  TYPE( nlls_control_type ) :: coptions
  type( C_PTR), value :: cw

  type( params_wrapper ) :: fparams
  TYPE( f_nlls_inform_type ) :: finform
  TYPE( f_nlls_workspace ), pointer :: fw
  TYPE( f_nlls_control_type ) :: foptions

  logical :: f_arrays

  ! copy data in and associate pointers correctly
  call copy_control_in(coptions, foptions, f_arrays)

  call c_f_procpointer(r, fparams%r)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  call c_f_pointer(cw, fw)
  fparams%params = params

  call f_ral_nlls_iterate( n, m, cx, &
       c_eval_r, c_eval_j,   &
       c_eval_hf, fparams,   &
       finform, foptions, fw)

  ! Copy data out
  call copy_info_out(finform, cinform)

end subroutine ral_nlls_iterate_d
