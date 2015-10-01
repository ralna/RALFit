program driver

use nlls_module

integer                   :: n, m, len_work_int, len_work_real
real(wp), allocatable     :: X(:), Work_real(:)
integer, allocatable      :: Work_int(:)
type( NLLS_inform_type )  :: status
type( NLLS_control_type ) :: options

write(*,*) '==============='
write(*,*) 'RAL NLLS driver'
write(*,*) '==============='
write(*,*) ' '


n = 100

m = 90

len_work_int = n
allocate( Work_int(len_work_int) )

len_work_real = n
allocate( Work_real(len_work_real) ) 

allocate( X(n) )

X = 1.0_wp

options%print_level = 3

call ral_nlls(n, m, X, Work_int, len_work_int, & 
              Work_real, len_work_real,         &
              eval_F, eval_J,                   &
              status, options )

end program driver


