! Copyright (C) 2016 Science and Technology Facilities Council (STFC).
! All rights reserved.
! examples/Fortran/Lanczos.f90 (based on)
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! examples/Fortran/LanczosFD.f90

module lanczos_module_fd
   use ral_nlls_double, only : params_base_type
   implicit none

   integer, parameter, public :: wp = kind(0d0)
   real(wp), parameter :: tol = 5.0e-2_wp

   type, extends(params_base_type) :: params_type
      real(wp), dimension(:), allocatable :: t ! The m data points t_i
      real(wp), dimension(:), allocatable :: y ! The m data points y_i
   end type params_type

contains

   subroutine eval_r(status, n, m, x, r, params)
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

   end subroutine eval_r

end module lanczos_module_fd


program lanczos_fd

   use ral_nlls_double, only: nlls_options, nlls_inform, ral_nlls_eval_j_dummy,&
      ral_nlls_eval_hf_dummy, nlls_solve
   use lanczos_module_fd

   implicit none

   type(nlls_options) :: options
   type(nlls_inform) :: inform

   integer :: m,n
   real(wp), allocatable :: x(:)
   type(params_type) :: params
   real(wp) :: tic, toc

   real(wp), parameter, Dimension(6) :: x_exp = (/ 8.6811579049232354E-002, &
      0.95495550838468568, 0.84399029667347680, 2.9515586935011613,         &
      1.5825909814491355, 4.9863421286471254 /)
   logical :: ok
   logical :: oki
   integer :: i

   Continue

   ! data to be fitted
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

   ! call fitting routine
   n = 6
   allocate(x(n))
   x = (/ 1.2, 0.3, 5.6, 5.5, 6.5, 7.6 /) ! SP 1
   ok = .False.

   options%print_level = 4
   options%use_ews_subproblem = .true.
   options%type_of_method = 2
   options%regularization_term = 1.0e-2
   options%regularization_power = 2.0
   options%reg_order = -1.0
   options%inner_method = 2
   options%maxit = 1000
   options%fd_step = 5.0e-3


   call cpu_time(tic)
   call nlls_solve(n,m,x,eval_r, ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, params, options, inform)
   if(inform%status.ne.0) then
      print *, "ral_nlls() returned with nonzero flag: ", inform%status
      goto 100
   endif
   call cpu_time(toc)

   ! Print result and check solution
   ok = .True.
   Write(*,*) 'Solution: '
   Write(*,Fmt=99998) 'idx',  'x',  'x*', 'error'
   Do i = 1, n
      oki = abs(x(i) - x_exp(i)) <= tol
      ok = ok .And. oki
      Write(*,Fmt=99999) i, x(i), x_exp(i), merge('PASS', 'FAIL', oki), abs(x(i) - x_exp(i))
   End Do
   print *, ""

   ! Print result
   print *, "Took ", inform%iter, " iterations"
   print *, "     ", inform%f_eval, " function evaluations"
   print *, "     ", inform%g_eval, " gradient evaluations (FD)"
   write(*,'(9X,Es9.2,2X,A)') toc-tic, "seconds"

100 Continue

   if (allocated(x)) deallocate(x)
   if (allocated(params%t)) deallocate(params%t)
   if (allocated(params%y)) deallocate(params%y)
   stop merge(0, 5, ok)

99999 Format (5X,I3,1X,2(Es13.6e2,2X),A4,2X,Es9.2e2)
99998 Format (5X,A3,1X,2(A13,2X),4X,2X,A9)
end program lanczos_fd
