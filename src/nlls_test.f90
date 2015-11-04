program nlls_test
  
! Test deck for nlls_module

  use nlls_module 
  implicit none

  real(wp), allocatable :: x(:),y(:),z(:)
  real(wp), allocatable :: A(:,:), B(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha
  integer :: m, n, i, no_errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test the helper subroutines !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  no_errors = 0

!! mult_J
  
  n = 2
  m = 4
  
  allocate(z(m*n),x(m),y(n))
  x = 1.0_wp
  z = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /)
  call mult_J(z,m,n,x,y)
  if ( norm2( y - (/16.0, 20.0 /) ) > 1e-12) then
     write(*,*) 'error :: mult_J test failed'
     no_errors = no_errors + 1 
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
     no_errors = no_errors + 1 
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
     no_errors = no_errors + 1     
  end if
  
  deallocate(x,A,B,results)
  
!! max_eig

  n = 4
  allocate(x(n),A(n,n), B(n,n))
  A = reshape( (/1.0, 2.0, 3.0, 4.0, &
                 2.0, 4.0, 6.0, 8.0, &
                 3.0, 6.0, 9.0, 12.0, & 
                 4.0, 8.0, 12.0, 16.0/), shape(A))
  B = 0.0_wp
  do i = 1,n
     B(i,i) = real(i,wp)
  end do

  call max_eig(A,B,n,alpha,x)

  if ( abs( alpha - 10.0 ) > 1e-12 ) then
     write(*,*) 'error :: max_eig test failed'
     no_errors = no_errors + 1 
  end if
  
  deallocate(A,B,x)

! Report back results....
  
  if (no_errors == 0) then
     write(*,*) '*** All tests passed successfully! ***'
  else
     write(*,*) 'There were ', no_errors,' errors'
  end if

end program nlls_test
