program nlls_test
  
! Test deck for nlls_module

  use nlls_module 
  implicit none

  real(wp), allocatable :: x(:)
  real(wp), allocatable :: A(:,:), B(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha
  integer :: n, i, no_errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test the helper subroutines !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  no_errors = 0

! outer_product
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
  
! max_eig

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

  write(*,*) 'alpha = ', alpha
  
! Report back results....
  
  if (no_errors == 0) then
     write(*,*) '*** All tests passed successfully! ***'
  else
     write(*,*) 'There were ', no_errors,' errors'
  end if

end program nlls_test
