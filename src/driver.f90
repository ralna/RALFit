program driver

use nlls_module
use example_module
implicit none

integer                   :: n, m, len_work_int, len_work_real, i
real(wp), allocatable     :: X(:), Work_real(:)
integer, allocatable      :: Work_int(:)
type( NLLS_inform_type )  :: status
type( NLLS_control_type ) :: options
type( user_type ), target :: params

write(*,*) '==============='
write(*,*) 'RAL NLLS driver'
write(*,*) '==============='
write(*,*) ' '

n = 2

m = 67

len_work_int = n
allocate( Work_int(len_work_int) )

len_work_real = n
allocate( Work_real(len_work_real) ) 

allocate( X(n) )

X(1) = 1.0 
X(2) = 2.0

options%print_level = 3

! Get params for the function evaluations
allocate(params%x_values(m))
allocate(params%y_values(m))

call generate_data_example(params%x_values,params%y_values,m)

call ral_nlls(n, m, X,                 &
              eval_F, eval_J, params,  &
              status, options )

do i = 1,n
   write(*,*) X(i)
end do

X(1) = 1.0 
X(2) = 2.0

call ral_nlls_int_func(n, m, X, Work_int, len_work_int, & 
              Work_real, len_work_real,         &
              status, options )

do i = 1,n
   write(*,*) X(i)
end do


end program driver

