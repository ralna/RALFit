! Copyright (C) 2016 Science and Technology Facilities Council (STFC).
! All rights reserved.
! examples/Fortran/Lanczos.f90 (based on)
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! examples/Fortran/LanczosFD2.f90

! This example showcases how to call the Finite Difference machinery
! explicitly instead of using the automatic implicit framework from the
! RALFit solver (i.e. passing dummy jacobian eval_j callback.)
! The example shows the usage of jacobian_setup, jacobian_calc and
! jacobian_free while solving the Lanczos problem. In this variant,
! the jacobian matrix is estimated every 10 calls to eval_j and for the
! next call only the exact jacobian is returned and this regime is repeated
! until conversion.

module lanczos_module_fd2

#if SINGLE_PRECISION
   use ral_nlls_single
#else
   use ral_nlls_double
#endif

   implicit none

   type, extends(params_base_type) :: params_type
      real(wp), dimension(:), allocatable :: t ! The m data points t_i
      real(wp), dimension(:), allocatable :: y ! The m data points y_i
      real(wp), dimension(:), allocatable :: r ! Auxiliary residual vector for FD
      type(jacobian_handle) :: handle ! The FD handle
      integer :: cnt, reset ! counter and reset
      real(wp) :: fd_step
      ! pointer to x(:)
      real(wp), dimension(:), pointer :: x
      ! telemetry
      integer :: f_cnt = 0, j_cnt = 0
   end type params_type

contains

   subroutine assign(p, t)
    Implicit None
    Real(Kind=wp), Dimension(:), target :: t
    Real(Kind=wp), Dimension(:), pointer :: p
    continue
    p => t
   end subroutine assign

   subroutine eval_r(status, n, m, x, r, params)
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

   end subroutine eval_r

   ! the call-back provides exact jacobian every CNT calls estimates t
   subroutine eval_J(status, n, m, x, J, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: J
      class(params_base_type), intent(inout) :: params

      continue

      status = 0

      select type(params)
       type is(params_type)
         params%cnt = params%cnt - 1
         if (params%cnt == 0) then
            params%cnt = params%reset
            J(    1:  m) = -exp(-x(2)*params%t(1:m))                     ! J_i1
            J(  m+1:2*m) = params%t(1:m) * x(1) * exp(-x(2)*params%t(1:m))! J_i2
            J(2*m+1:3*m) = -exp(-x(4)*params%t(1:m))                     ! J_i3
            J(3*m+1:4*m) = +params%t(1:m) * x(3) * exp(-x(4)*params%t(1:m))! J_i4
            J(4*m+1:5*m) = -exp(-x(6)*params%t(1:m))                     ! J_i5
            J(5*m+1:6*m) = +params%t(1:m) * x(5) * exp(-x(n)*params%t(1:m))! J_i6
            params%j_cnt = params%j_cnt + 1
         else
            ! get the residuals for the current iterate
            Call eval_r(status, n, m, x, params%r, params)
            params%f_cnt = params%f_cnt + 1
            if (status /= 0) return
            call assign(params%x, x(1:n))
            Call jacobian_calc(status, params%handle, params%x(1:n), params%r(1:m), &
               J(1:n*m), params%fd_step)
            params%f_cnt = params%f_cnt + n ! FD issues n calls to eval_r
         end if
      end select

   end subroutine eval_J

end module lanczos_module_fd2

program lanczos_fd2

   use lanczos_module_fd2

   implicit none

   type(nlls_options) :: options
   type(nlls_inform) :: inform

   integer :: m,n
   real(wp), allocatable, target :: x(:)

   type(params_type) :: params
   real(wp) :: tic, toc

   real(wp), parameter, Dimension(6) :: x_exp = (/ 8.6811579049232354E-002, &
      0.95495550838468568, 0.84399029667347680, 2.9515586935011613,         &
      1.5825909814491355, 4.9863421286471254 /)
   real(wp) :: tol = 0.0_wp
   logical :: ok
   logical :: oki
   integer :: i, status

   Continue

   ! data to be fitted
   m = 24
   allocate(params%t(m), params%y(m), params%r(m))
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
   ok = .False.

   options%print_level = 4
   options%use_ews_subproblem = .true.
   options%type_of_method = 2
   options%regularization_term = 1.0e-2
   options%regularization_power = 2.0
   options%reg_order = -1.0
   options%inner_method = 2
   options%maxit = 1000

   if (wp == lp) then
      ! Start solver closer to the expected solution when using low precision
      x(1:6) = (/ 8.68E-2, 0.95, 0.843, 2.95, 1.58, 4.98 /)
      options%fd_step = 1.0e-3
      tol = 5.0e-3_wp
      params%reset = 3 ! times FD is used before providing exact derivatives
   else
      x = (/ 1.2, 0.3, 5.6, 5.5, 6.5, 7.6 /) ! SP 1
      options%fd_step = 1.0e-5
      tol = 5.0e-4_wp
      params%reset = 10 ! times FD is used before providing exact derivatives
   end if

   ! setup up FD machinery

   params%cnt = params%reset
   params%fd_step = options%fd_step



   call jacobian_setup(status, params%handle, n, m, x, eval_r, params)
   if (status /= 0) then
      print *, 'Problems while calling jacobian_setup'
      goto 100
   endif

   call cpu_time(tic)
   call nlls_solve(n, m, x, eval_r, eval_j, ral_nlls_eval_hf_dummy, params, &
      options, inform)
   call jacobian_free(params%handle)
   if(inform%status.ne.0) then
      print *, "ral_nlls() returned with nonzero flag: ", inform%status
      goto 100
   endif
   call cpu_time(toc)

   ! Print result and check solution
   ok = .True.
   Write(*,*) 'Solution: '
   Write(*,Fmt=99998) 'idx',  'x',  'x*'
   Do i = 1, n
      oki = abs(x(i) - x_exp(i)) <= tol
      ok = ok .And. oki
      Write(*,Fmt=99999) i, x(i), x_exp(i), merge('PASS', 'FAIL', oki)
   End Do
   print *, ""

   ! Print result
   print *, "Took ", inform%iter, " iterations"
   print *, "     ", inform%f_eval + params%f_cnt, " (", params%f_cnt, ") function evaluations (FD)"
   print *, "     ", inform%g_eval + params%j_cnt, " (", params%j_cnt, ") gradient evaluations (FD)"
   write(*,'(9X,Es9.2,2X,A)') toc-tic, "seconds"

100 Continue

   if (allocated(x)) deallocate(x)
   if (allocated(params%t)) deallocate(params%t)
   if (allocated(params%y)) deallocate(params%y)
   if (allocated(params%r)) deallocate(params%r)

   stop merge(0, 5, ok)

99999 Format (5X,I3,1X,2(Es13.6e2,2X),A4)
99998 Format (5X,A3,1X,2(A13,2X))
end program lanczos_fd2
