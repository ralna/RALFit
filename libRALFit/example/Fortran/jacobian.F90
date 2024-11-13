! Copyright (c) 2017, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

! This example showcases how to call the Finite Difference machinery,
! illustrates the usage of jacobian_setup, jacobian_calc and
! jacobian_free as a stand-alone use case.

module jacobian_module_fd

#if SINGLE_PRECISION
   use ral_nlls_single
#else
   use ral_nlls_double
#endif

   implicit none

   type, extends(params_base_type) :: params_type
      real(wp), dimension(:), allocatable :: t ! The m data points t_i
      real(wp), dimension(:), allocatable :: y ! The m data points y_i
   end type params_type

contains

   subroutine eval_f(status, n, m, x, r, params)
      implicit none
      ! r_i = y_i - x_1 e^(-x_2 t_i) - x_3 e^(-x_4 t_i) - x_5 e^(-x_6 t_i)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: r
      class(params_base_type), intent(inout) :: params

      select type(params)
       type is(params_type)
         r(1:m) = params%y(:) &
            - x(1)*exp(-x(2)*params%t(:)) &
            - x(3)*exp(-x(4)*params%t(:)) &
            - x(5)*exp(-x(n)*params%t(:))
      end select

      status = 0 ! success

   end subroutine eval_f

   ! Exact Jacobian used to validate the FD Jacobian
   subroutine eval_J(status, n, m, x, J, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: J
      class(params_base_type), intent(inout) :: params

      select type(params)
       type is(params_type)
         J(    1:  m) = -exp(-x(2)*params%t(1:m))                     ! J_i1
         J(  m+1:2*m) = params%t(1:m) * x(1) * exp(-x(2)*params%t(1:m))! J_i2
         J(2*m+1:3*m) = -exp(-x(4)*params%t(1:m))                     ! J_i3
         J(3*m+1:4*m) = +params%t(1:m) * x(3) * exp(-x(4)*params%t(1:m))! J_i4
         J(4*m+1:5*m) = -exp(-x(6)*params%t(1:m))                     ! J_i5
         J(5*m+1:6*m) = +params%t(1:m) * x(5) * exp(-x(n)*params%t(1:m))! J_i6
      end select

      status = 0 ! Success
   end subroutine eval_J

end module jacobian_module_fd

program jacobian

   use jacobian_module_fd

   implicit none


   integer :: m, n
   real(wp), allocatable :: x(:), fx(:), j(:), j_exp(:)
   type(params_type) :: params
   real(wp) :: tic, toc, nrm2, tol, fd_step
   Type(jacobian_handle) :: handle
   logical :: ok, oki
   integer :: i, status, colj_start, colj_end

   Continue
   if (wp == lp) then
      ! Single precision
      tol = 7.0e-4_wp
      fd_step = 1.0e-3_wp
   else
      ! Double precision
      tol = 5.0e-7_wp
      fd_step = 1.0e-7_wp
   end if

   ! data to be fitted
   n = 6
   m = 24
   allocate(params%t(m), params%y(m))
   ! Data from Lanczos 3
   params%t(:) = (/ 0.00000E+00, &
      5.00000E-02, &
      1.00000E-01, &
      1.50000E-01, &
      2.00000E-01, &
      2.50000E-01, &
      3.00000E-01, &
      3.50000E-01, &
      4.00000E-01, &
      4.50000E-01, &
      5.00000E-01, &
      5.50000E-01, &
      6.00000E-01, &
      6.50000E-01, &
      7.00000E-01, &
      7.50000E-01, &
      8.00000E-01, &
      8.50000E-01, &
      9.00000E-01, &
      9.50000E-01, &
      1.00000E+00, &
      1.05000E+00, &
      1.10000E+00, &
      1.15000E+00 /)
   params%y(:) = (/ 2.5134E+00, &
      2.0443E+00, &
      1.6684E+00, &
      1.3664E+00, &
      1.1232E+00, &
      0.9269E+00, &
      0.7679E+00, &
      0.6389E+00, &
      0.5338E+00, &
      0.4479E+00, &
      0.3776E+00, &
      0.3197E+00, &
      0.2720E+00, &
      0.2325E+00, &
      0.1997E+00, &
      0.1723E+00, &
      0.1493E+00, &
      0.1301E+00, &
      0.1138E+00, &
      0.1000E+00, &
      0.0883E+00, &
      0.0783E+00, &
      0.0698E+00, &
      0.0624E+00 /)

   allocate(x(n), fx(m), j(n*m), j_exp(n*m))
   x = (/ 1.2, 0.3, 5.6, 5.5, 6.5, 7.6 /) ! SP 1
   x(:) = 1.0_wp

   call cpu_time(tic)

   ! setup up FD machinery and set storace scheme to fortran (column-major)
   call jacobian_setup(status, handle, n, m, x, eval_f, params, f_storage=.True.)
   if (status /= 0) then
      print *, 'Problem while calling jacobian_setup'
      goto 100
   endif

   ! evaluate the residuals at point x
   call eval_f(status, n, m, x, fx, params)
   if (status /= 0) then
      print *, 'Problem while calling eval_f'
      goto 100
   endif

   ! estimate the jacobian matrix at point x, fx is the residual at point x
   call jacobian_calc(status, handle, x, fx, j, fd_step=fd_step)
   if (status /= 0) then
      print *, 'Problem while calling jacobian_calc'
      goto 100
   endif

   call jacobian_free(handle)
   call cpu_time(toc)

   ! Print result and check solution
   call eval_j(status, n, m, x, j_exp, params)
   if (status /= 0) then
      print *, 'Problem while calling eval_j'
      goto 100
   endif

   ok = .True.
   print *, ""
   Write(*,Fmt=99998) 'col-idx',  'norm(Ji-Ji*)'
   Do i = 1, n
      colj_start = (i-1)*m + 1
      colj_end = i * m
      nrm2 = norm2(j(colj_start:colj_end) - j_exp(colj_start:colj_end))
      oki = nrm2 <= tol
      ok = ok .And. oki
      Write(*,Fmt=99999) i, nrm2, merge('PASS', 'FAIL', oki)
   End Do
   print *, ""

   ! Print result
   print *, "Took "
   write(*,'(9X,Es9.2,2X,A)') toc-tic, "seconds"

100 Continue

   if (allocated(x)) deallocate(x)
   if (allocated(fx)) deallocate(fx)
   if (allocated(j)) deallocate(j)
   if (allocated(j_exp)) deallocate(j_exp)
   if (allocated(params%t)) deallocate(params%t)
   if (allocated(params%y)) deallocate(params%y)

   stop merge(0, 5, ok)

99999 Format (5X,I3,1X,1(Es13.6e2,2X),A4)
99998 Format (5X,A3,1X,1(A13,2X))
end program jacobian
