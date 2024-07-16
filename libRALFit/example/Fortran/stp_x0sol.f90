! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! Copyright (c) 2019, The Numerical Algorithms Group Ltd (NAG)
! All rights reserved.
! Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! examples/Fortran/Lanczos.f90
! STP to test that initial point is solution.

module lanczos_module
  use ral_nlls_double, only : params_base_type
  implicit none

  integer, parameter :: wp = kind(0d0)

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
            - x(5)*exp(-x(6)*params%t(:))
    end select

    status = 0 ! success

  end subroutine eval_r

  subroutine eval_J(status, n, m, x, J, params)
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
         J(5*m+1:6*m) = +params%t(1:m) * x(5) * exp(-x(6)*params%t(1:m))! J_i6
    end select


    status = 0 ! Success
  end subroutine eval_J

  subroutine eval_HF(status, n, m, x, r, HF, params)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(in) :: r
      real(wp), dimension(*), intent(out) :: HF
      class(params_base_type), intent(inout) :: params

      HF(1:n*n) = 0.0
      select type(params)
      type is(params_type)
         HF(    2) = sum( r(1:m) * params%t(1:m) * exp(-x(2)*params%t(1:m)))  ! H_21
         HF(1*n+1) = HF(2)                                                     ! H_12
         HF(1*n+2) = sum(-r(1:m) * (params%t(1:m)**2) * x(1) * exp(-x(2)*params%t(1:m)))! H_22
         HF(2*n+4) = sum( r(1:m) * params%t(1:m) * exp(-x(4)*params%t(1:m)))  ! H_43
         HF(3*n+3) = HF(2*n+4)                                                 ! H_34
         HF(3*n+4) = sum(-r(1:m) * (params%t(1:m)**2) * x(3) * exp(-x(4)*params%t(1:m)))! H_44
         HF(4*n+6) = sum( r(1:m) * params%t(1:m) * exp(-x(6)*params%t(1:m)))  ! H_65
         HF(5*n+5) = HF(4*n + 6)                                                     ! H_56
         HF(5*n+6) = sum(-r(1:m) * (params%t(1:m)**2) * x(5) * exp(-x(6)*params%t(1:m)))! H_66
      end select

      status = 0 ! Success
    end subroutine eval_HF

    subroutine eval_HP(status, n, m, x, y, HP, params)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(in) :: y
      real(wp), dimension(*), intent(out) :: HP
      class(params_base_type), intent(inout) :: params

      integer :: i

      HP(1:n*m) = 0.0
      select type(params)
      type is (params_type)
         do i = 1, m
            HP((n-1)*i + 1) = params%t(i) * exp(-x(2)*params%t(i)) * y(2)
            HP((n-1)*i + 2) = params%t(i) * exp(-x(2)*params%t(i)) * y(1) - &
                 (params%t(i)**2) * x(1) * exp(-x(2)*params%t(i)) * y(2)
            HP((n-1)*i + 3) = params%t(i) * exp(-x(4)*params%t(i)) * y(4)
            HP((n-1)*i + 4) = params%t(i) * exp(-x(4)*params%t(i)) * y(3) - &
                 (params%t(i)**2) * x(3) * exp(-x(4)*params%t(i)) * y(4)
            HP((n-1)*i + 5) = params%t(i) * exp(-x(6)*params%t(i)) * y(6)
            HP((n-1)*i + 6) = params%t(i) * exp(-x(6)*params%t(i)) * y(5) - &
                 (params%t(i)**2) * x(5) * exp(-x(6)*params%t(i)) * y(6)
         end do
      end select


    end subroutine eval_HP


end module lanczos_module


program lanczos

  use ral_nlls_double
  use lanczos_module

  implicit none

  type(nlls_options) :: options
  type(nlls_inform) :: inform

  integer :: m,n
  real(wp), allocatable :: x(:)
  type(params_type) :: params
  integer :: inner_method
  real(wp) :: tic, toc
  logical :: ok
  continue
  ok = .False.
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
  x=(/ 8.6810580646277377E-02_wp, 0.9549498952571175_wp, 0.8439878897334759_wp,&
       2.9515524061595775_wp, 1.5825943863392782_wp, 4.9863401266962830_wp /)

  options%print_level = 1
  options%exact_second_derivatives = .true.
  options%model = 4
  options%nlls_method = 3
  options%use_ews_subproblem = .true.
  options%type_of_method = 2
  options%regularization_term = 1.0e-2
  options%regularization_power = 2.0
  options%reg_order = -1.0
  options%inner_method = 2
  options%maxit = 10

  call cpu_time(tic)
  call nlls_solve(n,m,x,eval_r, eval_J, eval_HF, params, options, inform)!, eval_HP=eval_HP)
  call cpu_time(toc)
  if(inform%status.ne.0) then
     print *, "nlls_solve returned with error flag: ", inform%status
     goto 100
  endif

  ! Print result
  print *, "Found a local optimum at x = ", x
  print *, "Took ", inform%iter, " iterations"
  print *, "     ", inform%f_eval, " function evaluations"
  print *, "     ", inform%g_eval, " gradient evaluations"
  print *, "     ", inform%h_eval, " hessian evaluations"
  print *, "     ", toc-tic, " seconds"
  ok = inform%iter == 0 .And. inform%f_eval == 1 .And. inform%g_eval == 1

100 Continue

  If (Allocated(params%t)) Deallocate(params%t)
  If (Allocated(params%y)) Deallocate(params%y)

  Stop merge(0, 7, ok)
end program lanczos
