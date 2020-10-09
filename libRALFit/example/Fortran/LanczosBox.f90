! examples/Fortran/LanczosBox.f90

module lanczos_box_module
  use ral_nlls_double
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
      status = 0 ! success
      return
    Class Default
      stop 'evalr: ERROR do not know how to handle this class...'
    end select
    status = -1 ! fail
  end subroutine eval_r

  subroutine eval_J(status, n, m, x, J, params)
    integer, intent(out) :: status
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(wp), dimension(*), intent(in) :: x
    real(wp), dimension(*), intent(out) :: J
  class(params_base_type), intent(inout) :: params
    Integer :: r, c

    select type(params)
    type is(params_type)
      J(    1:  m) = -exp(-x(2)*params%t(1:m))                       ! J_i1
      J(  m+1:2*m) = params%t(1:m) * x(1) * exp(-x(2)*params%t(1:m)) ! J_i2
      J(2*m+1:3*m) = -exp(-x(4)*params%t(1:m))                       ! J_i3
      J(3*m+1:4*m) = +params%t(1:m) * x(3) * exp(-x(4)*params%t(1:m))! J_i4
      J(4*m+1:5*m) = -exp(-x(6)*params%t(1:m))                       ! J_i5
      J(5*m+1:6*m) = +params%t(1:m) * x(5) * exp(-x(6)*params%t(1:m))! J_i6
      status = 0 ! success
      return
    Class Default
      stop 'evalr: ERROR do not know how to handle this class...'
    end select
    status = -1 ! fail
  end subroutine eval_J

  subroutine eval_HF(status, n, m, x, r, HF, params)
    integer, intent(out) :: status
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(wp), dimension(*), intent(in) :: x
    real(wp), dimension(*), intent(in) :: r
    real(wp), dimension(*), intent(out) :: HF
  class(params_base_type), intent(inout) :: params

    HF(1:n*n) = 0
    select type(params)
    type is(params_type)
      HF(    2) = sum( r(1:m) * params%t(1:m) * exp(-x(2)*params%t(1:m)))            ! H_21
      HF(1*n+1) = HF(2)                                                              ! H_12
      HF(1*n+2) = sum(-r(1:m) * (params%t(1:m)**2) * x(1) * exp(-x(2)*params%t(1:m)))! H_22
      HF(2*n+4) = sum( r(1:m) * params%t(1:m) * exp(-x(4)*params%t(1:m)))            ! H_43
      HF(3*n+3) = HF(2*n+4)                                                          ! H_34
      HF(3*n+4) = sum(-r(1:m) * (params%t(1:m)**2) * x(3) * exp(-x(4)*params%t(1:m)))! H_44
      HF(4*n+6) = sum( r(1:m) * params%t(1:m) * exp(-x(6)*params%t(1:m)))            ! H_65
      HF(5*n+5) = HF(4*n + 6)                                                        ! H_56
      HF(5*n+6) = sum(-r(1:m) * (params%t(1:m)**2) * x(5) * exp(-x(6)*params%t(1:m)))! H_66
      status = 0 ! success
      return
    Class Default
      stop 'evalr: ERROR do not know how to handle this class...'
    end select
    status = -1 ! fail
  end subroutine eval_HF
  
end module lanczos_box_module

program lanczos_box
  use lanczos_box_module
  implicit none


  type(nlls_options) :: options
  type(nlls_inform) :: inform

  integer :: m,n
  real(wp), allocatable :: x(:), xnew(:), blx(:), bux(:)
  type(params_type) :: params
  integer :: i
  real(wp) :: tic, toc
 
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

  n = 6
  allocate(x(n))
  x = (/ 0.44, -0.2408, 2.598, 3.44, -6.199, 0.0977 /)
 
  ! Add bounds on the variables
  Allocate(blx(n),bux(n),xnew(n))
  blx(1:n) = -1.0
  bux(1:n) = +1.0e+20
  blx(1) = 0.445
  bux(1) = 1.0
  blx(5) = -1.0
  bux(5) = 1.0
  bux(6) = 10.0

  options%print_level = 3
  options%maxit = 100
  options%exact_second_derivatives = .true.
  options%model = 3
  options%type_of_method = 1 ! TR / Reg
  options%inner_method = 1 ! Implicit / min / basereg
  options%nlls_method = 4 ! Dogleg / AINT / More-Sorensen (LinSolve) / Galahad
  options%box_max_ntrfail = 2

!!$  Call nlls_setup_bounds(params, n, blx, bux, options, inform)
!!$  if (inform%status/=0) then
!!$    Write(*,*) 'ERROR: nlls_setup_bounds failed, status=', inform%status
!!$    stop
!!$  End if


  ! call fitting routine
  call cpu_time(tic)
  call nlls_solve(n,m,x,eval_r, eval_J, eval_HF, params, options, inform, &
       lower_bounds=blx, upper_bounds=bux)
  if(inform%status.ne.0) then
     print *, "ral_nlls() returned with error flag ", inform%status
     stop
  endif
  call cpu_time(toc)

  ! Print result
  Write(*,*) 'Solution: '
  Write(*,Fmt=99998) 'idx', 'low bnd', 'x', 'upp bnd'
  Do i = 1, n
     Write(*,Fmt=99999) i, blx(i), x(i), bux(i)
  End Do
  print *, ""
  print *, "Objective Value at solution    = ", inform%obj
  print *, "Objective Gradient at solution = ", inform%norm_g
  print *, "Took ", inform%iter, " iterations (LS:",inform%ls_step_iter,' PG:',inform%pg_step_iter,')'
  print *, "     ", inform%f_eval, " function evaluations (LS:",inform%f_eval_ls,' PG:',inform%f_eval_pg,')'
  print *, "     ", inform%g_eval, " gradient evaluations (LS:",inform%g_eval_ls,' PG:',inform%g_eval_pg,')'
  print *, "     ", inform%h_eval, " hessian (eval Hf) evaluations"
  print *, "     ", inform%hp_eval, " hessian (eval_HP) evaluations"
  print *, "     ", toc-tic, " seconds"

99999  Format (5X,I3,1X,3(Es13.6e2,2X))
99998  Format (5X,A3,1X,3(A13,2X))
99997  Format (5X,I3,16X,Es13.6e2)
end program lanczos_box
