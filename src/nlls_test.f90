program nlls_test
  
! Test deck for nlls_module

  use nlls_module 
  use example_module
  implicit none

  type( NLLS_inform_type )  :: status
  type( NLLS_control_type ) :: options
  type( user_type ), target :: params
  real(wp), allocatable :: x(:),y(:),z(:)
  real(wp), allocatable :: A(:,:), B(:,:), C(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha
  integer :: m, n, i, no_errors_helpers, no_errors_main, info
  integer :: nlls_method

  type( NLLS_workspace ) :: work, work2

!!!!!!!!!!!!!!!!!!!!!!!!
!! Test the main file !!
!!!!!!!!!!!!!!!!!!!!!!!!

no_errors_main = 0

n = 2

m = 67

do nlls_method = 1,2
   allocate( x(n) )

   X(1) = 1.0 
   X(2) = 2.0
   
   options%print_level = 0
   options%nlls_method = nlls_method

   ! Get params for the function evaluations
   allocate(params%x_values(m))
   allocate(params%y_values(m))
   
   call generate_data_example(params%x_values,params%y_values,m)
   
   call ral_nlls(n, m, X,                         &
        eval_F, eval_J, eval_H, params,  &
        status, options )
   if ( status%status .ne. 0 ) then
      write(*,*) 'ral_nlls failed to converge:'
      write(*,*) 'NLLS_METHOD = ', nlls_method
      no_errors_main = no_errors_main + 1
   end if
   
   deallocate(x,params%x_values,params%y_values)
   
end do

if (no_errors_main == 0) then
   write(*,*) '*** All (main) tests passed successfully! ***'
else
   write(*,*) 'There were ', no_errors_main,' errors'
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test the helper subroutines !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
no_errors_helpers = 0

!! mult_J
  
  n = 2
  m = 4
  
  allocate(z(m*n),x(m),y(n))
  x = 1.0_wp
  z = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /)
  call mult_J(z,m,n,x,y)
  if ( norm2( y - (/16.0, 20.0 /) ) > 1e-12) then
     write(*,*) 'error :: mult_J test failed'
     no_errors_helpers = no_errors_helpers + 1 
  end if


  deallocate(z, x, y)

!! mult_Jt
  
  n = 2
  m = 4
  
  allocate(z(m*n),x(m),y(n))
  x = 1.0_wp
  z = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /)
  call mult_Jt(z,n,m,x,y)
  if ( norm2( y - (/10.0, 26.0 /) ) > 1e-12) then
     write(*,*) 'error :: mult_Jt test failed'
     no_errors_helpers = no_errors_helpers + 1 
  end if

  deallocate(z, x, y)

!! outer_product
  n = 4
  allocate(x(n), A(n,n), B(n,n), results(n))
  x = (/ 1.0, 2.0, 3.0, 4.0 /)
  A = reshape( (/1.0, 2.0, 3.0, 4.0, &
                 2.0, 4.0, 6.0, 8.0, &
                 3.0, 6.0, 9.0, 12.0, & 
                 4.0, 8.0, 12.0, 16.0/), shape(A))
  call outer_product(x,n,B)
  do i = 1, n
      results(i) = norm2(A(i,:) - B(i,:))
  end do
  if (norm2(results) > 1e-12) then
     write(*,*) 'error :: outer_product test failed'
     no_errors_helpers = no_errors_helpers + 1     
  end if
  
  deallocate(x,A,B,results)
  
!! max_eig
  
  n = 4
  m = 4
  ! make sure max_eig gets called
  options%nlls_method = 2
  call setup_workspaces(work,2,2,options,info) 

  allocate(x(n),A(n,n), B(n,n))
  A = reshape( (/1.0, 2.0, 3.0, 4.0, &
                 2.0, 4.0, 6.0, 8.0, &
                 3.0, 6.0, 9.0, 12.0, & 
                 4.0, 8.0, 12.0, 16.0/), shape(A))
  B = 0.0_wp
  do i = 1,n
     B(i,i) = real(i,wp)
  end do

  call max_eig(A,B,n,alpha,x,info,C,work%calculate_step_ws%AINT_tr_ws%max_eig_ws)

  if ( (abs( alpha - 10.0 ) > 1e-12).or.(info .ne. 0) ) then
     write(*,*) 'error :: max_eig test failed'
     no_errors_helpers = no_errors_helpers + 1 
  end if
  
  deallocate(A,B,x)

  ! check the 'hard' case...
  n = 4
  allocate(x(n),A(n,n), B(n,n))
  A = 0.0_wp  
  A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
  A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
  B = A
  A(1,1) = 1.0_wp; A(2,2) = 1.0_wp
  
  call max_eig(A,B,n,alpha,x,info,C,work%calculate_step_ws%AINT_tr_ws%max_eig_ws)

  if (.not. allocated(C)) then ! check C returned 
     write(*,*) 'error :: hard case of max_eig test failed - C not returned'
     no_errors_helpers = no_errors_helpers + 1 
  else
     allocate(y(2))
     y = shape(C)
     if ((y(1) .ne. 2) .or. (y(2) .ne. n)) then
        write(*,*) 'error :: hard case of max_eig test failed - wrong shape C returned'
        no_errors_helpers = no_errors_helpers + 1 
     else
        allocate(results(n))
        ! Repopulate A (was overwritten by eig routine)
        A = 0.0_wp  
        A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
        A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
        B = A
        A(1,1) = 1.0_wp; A(2,2) = 1.0_wp
        do i = 1, n
           results(i) = norm2(                        &
                matmul( A(3:4,3:4),C(1:2,i) )         &
                - alpha * matmul(B(3:4,3:4),C(1:2,i)) & 
                )
        end do
        if (norm2(results) > 1e-10) then
           write(*,*) 'error :: hard case of max_eig test failed - wrong vectors returned'
           write(*,*) 'results = ', results
           no_errors_helpers = no_errors_helpers + 1 
        end if
     end if
  end if
  
  deallocate(A,B,C,x,y,results)
  

  call setup_workspaces(work2,1,1,options,info)  !todo: deallocation routine
  ! check the error return
  n = 2
  allocate(x(n), A(n,n), B(n,n))
  A = 0.0_wp
  B = 0.0_wp
  A(1,2) = 1.0_wp
  A(2,1) = -1.0_wp
  B(1,1) = 1.0_wp
  B(2,2) = 1.0_wp

  call max_eig(A,B,n,alpha,x,info,C,work2%calculate_step_ws%AINT_tr_ws%max_eig_ws)
  if (info .ne. 1) then
     write(*,*) 'error :: all complex part of max_eig test failed'
     no_errors_helpers = no_errors_helpers + 1
  end if
  
  call max_eig(A,B,n+1,alpha,x,info,C,work2%calculate_step_ws%AINT_tr_ws%max_eig_ws)
  if (info .ne. 2) then
     write(*,*) 'error :: even part of max_eig test failed'
     no_errors_helpers = no_errors_helpers + 1
  end if
  
  deallocate(A,B,x)

!! matmult_outer
  
  n = 2
  m = 3
  allocate(A(m,n),B(m,m),C(m,m),results(m))
  A = reshape( (/1.0, 2.0, 3.0,  &
                 2.0, 4.0, 6.0/),&
                 shape(A))
  call matmult_outer(A,n,m,B)
  C = reshape( (/ 5.0, 10.0, 15.0,  &
       10.0, 20.0, 30.0, & 
       15.0, 30.0, 45.0 /) &
       , shape(C))
  do i = 1,m
     results(i) = norm2(C(:,i) - B(:,i))
  end do
  if (norm2(results) > 1e-10) then
     write(*,*) 'error :: matmult_outer test failed'
     no_errors_helpers = no_errors_helpers + 1
  end if
  
  deallocate(A,B,C,results)

 !! matmult_inner
  
  n = 2
  m = 3
  allocate(A(m,n),B(n,n),C(n,n),results(n))
  A = reshape( (/1.0, 2.0, 3.0,  &
                 2.0, 4.0, 6.0/),&
                 shape(A))
  call matmult_inner(A,n,m,B)
  C = reshape( (/ 14.0, 28.0,  &
                  28.0, 56.0 /) &
                  , shape(C))
  do i = 1,n
     results(i) = norm2(C(:,i) - B(:,i))
  end do
  if (norm2(results) > 1e-10) then
     write(*,*) 'error :: matmult_inner test failed'
     no_errors_helpers = no_errors_helpers + 1
  end if

  deallocate(A,B,C,results)
  

  ! test findbeta

  n = 3
  allocate(x(n),y(n),z(n))

  x = (/ 1.0, 2.0, 3.0 /) 
  y = (/ 2.0, 1.0, 1.0 /)
  
  call findbeta(x,y,1.0_wp,10.0_wp,alpha,info)

  if (info .ne. 0) then
     write(*,*) 'error -- findbeta did not work: info /= 0'
     no_errors_helpers = no_errors_helpers + 1
  else if ( ( norm2( x + alpha * y ) - 10.0_wp ) > 1e-12 ) then
     write(*,*) 'error -- findbeta did not work'
     write(*,*) '|| x + beta y|| = ', norm2( (x + alpha * y)-10.0_wp)
     no_errors_helpers = no_errors_helpers + 1
  end if
  
  deallocate(x,y,z)
  
! Report back results....
  
  if (no_errors_helpers == 0) then
     write(*,*) '*** All (helper) tests passed successfully! ***'
  else
     write(*,*) 'There were ', no_errors_helpers,' errors'
  end if


 if (no_errors_helpers == 0) then
    write(*,*) ' '
    write(*,*) '**************************************'
    write(*,*) '*** All tests passed successfully! ***'
    write(*,*) '**************************************'
  end if
end program nlls_test
