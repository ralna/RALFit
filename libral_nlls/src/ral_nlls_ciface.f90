module ral_nlls_ciface

  use iso_c_binding
  use nlls_module, only:          & 
       f_nlls_control_type => nlls_control_type, &
       f_nlls_inform_type  => nlls_inform_type, &
       f_ral_nlls => ral_nlls, & 
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
     real(wp) :: initial_radius
     real(wp) :: maximum_radius
     real(wp) :: eta_successful
     real(wp) :: eta_very_successful
     real(wp) :: eta_too_successful
     real(wp) :: radius_increase
     real(wp) :: radius_reduce
     real(wp) :: radius_reduce_max
  end type nlls_control_type

  type, bind(C) :: nlls_inform_type 
     integer(C_INT) :: status
  end type nlls_inform_type

  abstract interface
     subroutine c_eval_f_type(fstatus, n, m, x, f, params)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       integer(c_int) :: fstatus
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: f
       type(C_PTR), value :: params
     end subroutine c_eval_f_type
  end interface

  abstract interface
     subroutine c_eval_j_type(fstatus, n, m, x, j, params)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       integer(c_int) :: fstatus
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(out) :: j
       type(C_PTR), value :: params
     end subroutine c_eval_j_type
  end interface

  abstract interface
     subroutine c_eval_hf_type(fstatus, n, m, x, f, hf, params)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n,m
       integer(c_int) :: fstatus
       real(c_double), dimension(*), intent(in) :: x
       real(c_double), dimension(*), intent(in) :: f
       real(c_double), dimension(*), intent(out) :: hf
       type(C_PTR), value :: params
     end subroutine c_eval_hf_type
  end interface

  type, extends(f_params_base_type) :: params_wrapper
     procedure(c_eval_f_type), nopass, pointer :: f
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
    fcontrol%initial_radius = ccontrol%initial_radius
    fcontrol%maximum_radius = ccontrol%maximum_radius
    fcontrol%eta_successful = ccontrol%eta_successful
    fcontrol%eta_very_successful = ccontrol%eta_very_successful
    fcontrol%eta_too_successful = ccontrol%eta_too_successful
    fcontrol%radius_increase = ccontrol%radius_increase
    fcontrol%radius_reduce = ccontrol%radius_reduce
    fcontrol%radius_reduce_max = ccontrol%radius_reduce_max

  end subroutine copy_control_in

  subroutine copy_info_out(finfo,cinfo)

    type(f_nlls_inform_type), intent(in) :: finfo
    type(nlls_inform_type) , intent(out) :: cinfo
    
    cinfo%status = finfo%status
    
  end subroutine copy_info_out

  subroutine c_eval_f(evalfstatus, n, m, x, f, fparams)
    integer, intent(in) :: n, m 
    integer, intent(out) :: evalfstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(out) :: f  
    class(f_params_base_type), intent(in) :: fparams

    select type(fparams)
    type is(params_wrapper)
       call fparams%f(evalfstatus,n,m,x(1:n),f(1:m),fparams%params)
    end select
    
  end subroutine c_eval_f

  subroutine c_eval_j(evaljstatus, n, m, x, j, fparams)
    integer, intent(in) :: n, m 
    integer, intent(out) :: evaljstatus
    double precision, dimension(*), intent(in) :: x
    double precision, dimension(*), intent(out) :: j
    class(f_params_base_type), intent(in) :: fparams

    select type(fparams)
    type is(params_wrapper)
       call fparams%j(evaljstatus,n,m,x(1:n),j(1:n*m),fparams%params)
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
       call fparams%hf(evalhstatus,n,m,x(1:n),f(1:m),hf(1:n**2),fparams%params)
    end select
    
  end subroutine c_eval_hf
  
end module ral_nlls_ciface

subroutine nlls_default_control_d(ccontrol) bind(C)
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
  ccontrol%initial_radius = fcontrol%initial_radius
  ccontrol%maximum_radius = fcontrol%maximum_radius
  ccontrol%eta_successful = fcontrol%eta_successful
  ccontrol%eta_very_successful = fcontrol%eta_very_successful
  ccontrol%eta_too_successful = fcontrol%eta_too_successful
  ccontrol%radius_increase = fcontrol%radius_increase
  ccontrol%radius_reduce = fcontrol%radius_reduce
  ccontrol%radius_reduce_max = fcontrol%radius_reduce_max

  
end subroutine nlls_default_control_d

subroutine ral_nlls_d(n, m, cx,  &
                      f, j, hf,  & 
                      params,            &
                      cinform, coptions) &
                      bind(C)

  use ral_nlls_ciface 
  implicit none

  integer, INTENT( IN ), value :: n, m
  real( wp ), dimension(*) :: cx
  type( C_FUNPTR ), value :: f
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

  call c_f_procpointer(f, fparams%f)
  call c_f_procpointer(j, fparams%j)
  call c_f_procpointer(hf, fparams%hf)
  fparams%params = params

  call f_ral_nlls( n, m, cx, &
       c_eval_f, c_eval_j,   &
       c_eval_hf, fparams,   &
       finform, foptions)

  ! Copy data out
   call copy_info_out(finform, cinform)
  
end subroutine ral_nlls_d

