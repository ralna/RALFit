module ral_nlls_ciface

  use iso_c_binding
  use nlls_module, only:          & 
       f_nlls_control_type => nlls_control_type, &
       f_nlls_inform_type  => nlls_inform_type, &
       f_ral_nlls => ral_nlls, & 
       f_ral_nlls_int_func => ral_nlls_int_func
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
!!$
!!$subroutine ral_nlls_d(n, m, cX, &
!!$                      cWork_int,len_work_int, &
!!$                      cWork_real,len_work_real, &
!!$                      ceval_F, ceval_J, & 
!!$                      cstatus, coptions) &
!!$                      bind(C)
!!$
!!$  use ral_nlls_ciface 
!!$  implicit none
!!$
!!$  integer( C_INT ), INTENT( IN ) :: n, m, len_work_int, len_work_real
!!$  type( C_PTR ), value :: cX
!!$  type( C_PTR ), value :: cWork_int, cWork_real
!!$  TYPE( f_nlls_inform_type )  :: cstatus
!!$  TYPE( f_nlls_control_type ) :: coptions
!!$  external :: ceval_F, ceval_J
!!$
!!$  INTERFACE
!!$       SUBROUTINE eval_F( status, X, f )
!!$         USE ISO_FORTRAN_ENV
!!$         
!!$         INTEGER ( int32 ), INTENT( OUT ) :: status
!!$         REAL ( real64 ), DIMENSION( : ),INTENT( OUT ) :: f
!!$         REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X
!!$         
!!$       END SUBROUTINE eval_F
!!$    END INTERFACE
!!$
!!$    INTERFACE
!!$       SUBROUTINE eval_J( status, X, J )
!!$         USE ISO_FORTRAN_ENV
!!$
!!$         INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
!!$         INTEGER ( int32 ), INTENT( OUT ) :: status
!!$         REAL ( real64 ), DIMENSION( : ),INTENT( IN ) :: X
!!$         REAL ( real64 ), DIMENSION( : , : ),INTENT( OUT ) :: J
!!$       END SUBROUTINE eval_J
!!$    END INTERFACE
!!$  
!!$
!!$  REAL( wp ), DIMENSION( n ) :: fX
!!$  INTEGER( C_INT ) :: fWork_int(len_work_int)
!!$  REAL( wp )       :: fWork_real(len_work_real)
!!$  TYPE( nlls_inform_type ) :: fstatus
!!$  TYPE( nlls_control_type ) :: foptions
!!$
!!$  logical :: f_arrays
!!$
!!$  ! copy data in and associate pointers correctly
!!$  call copy_control_in(coptions, foptions, f_arrays)
!!$
!!$  call f_ral_nlls( n, m, fX, &
!!$       fWork_int,len_work_int, &
!!$       fWork_real,len_work_real, &
!!$       feval_F, feval_J, & 
!!$       fstatus, foptions)
!!$
!!$  ! Copy data out
!!$   call copy_info_out(fstatus, cstatus)
!!$  
!!$end subroutine ral_nlls_d


subroutine ral_nlls_int_func_d(n, m, cX, &
     cstatus, coptions) &
     bind(C)

  use ral_nlls_ciface 
  implicit none

  integer( C_INT ), value, INTENT( IN ) :: n, m
  type( C_PTR ), value :: cX
  TYPE( nlls_inform_type )  :: cstatus
  TYPE( nlls_control_type ) :: coptions

  REAL( wp ), DIMENSION( : ), pointer :: fX
  TYPE( f_nlls_inform_type ) :: fstatus
  TYPE( f_nlls_control_type ) :: foptions

  logical :: f_arrays
  
  ! copy data in and associate pointers correctly
  call copy_control_in(coptions, foptions, f_arrays)
  
  call C_F_POINTER(cX, fX, shape = (/ n /) )

  call f_ral_nlls_int_func( n, m, fX, &
       fstatus, foptions)

  ! Copy data out
   call copy_info_out(fstatus, cstatus)
  
end subroutine ral_nlls_int_func_d

