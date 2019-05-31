! examples/Fortran/LanczosBox.f90

module lanczos_box_module
  use ral_nlls_double
  implicit none
  integer, parameter :: wp = kind(0d0)

  type, extends(params_box_type) :: params_type
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
    
!Write(*,*) 'In eval_r, x=',x(1:n)
    select type(params)
    type is(params_type)
       r(1:m) = params%y(:) &
            - x(1)*exp(-x(2)*params%t(:)) &
            - x(3)*exp(-x(4)*params%t(:)) &
            - x(5)*exp(-x(6)*params%t(:)) 
      status = 0 ! success
!Write(*,*) 'Out eval_r, r=',r(1:m)
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
         J(    1:  m) = -exp(-x(2)*params%t(1:m))                     ! J_i1
         J(  m+1:2*m) = params%t(1:m) * x(1) * exp(-x(2)*params%t(1:m))! J_i2
         J(2*m+1:3*m) = -exp(-x(4)*params%t(1:m))                     ! J_i3
         J(3*m+1:4*m) = +params%t(1:m) * x(3) * exp(-x(4)*params%t(1:m))! J_i4
         J(4*m+1:5*m) = -exp(-x(6)*params%t(1:m))                     ! J_i5
         J(5*m+1:6*m) = +params%t(1:m) * x(5) * exp(-x(6)*params%t(1:m))! J_i6

!          Do r = 1, m
!           Do c = 1, n
!             Write(*,'(Es9.2e2,1X)',advance='no') J(n*(r-1)+c)
!           End Do
!           Write(*,*) ''
!          End Do

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

      HF(1:n*n) = 0
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
    
  
end module lanczos_box_module


program lanczos_box

  use lanczos_box_module
  
  implicit none

  type(nlls_options) :: options
  type(nlls_inform) :: inform

  integer :: m,n
  real(wp), allocatable :: x(:), xnew(:), d(:), blx(:), bux(:)
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

  ! call fitting routine
  n = 6
  allocate(x(n))
  x = (/ 1.2, 0.3, 5.6, 5.5, 6.5, 7.6 /) ! SP 1
  x = (/ 1.2, 0.3, 5.6, 5.5, 6.5, 7.6 /) ! SP 2
 
  ! Add the box bound
  Allocate(blx(n),bux(n),xnew(n), d(n))
! TEST ONE  
  blx(1:n) = -1.0e20_wp
  bux(1:n) = +1.0e20_wp
!   blx(1:n) = 0.0e20_wp
!   bux(1:n) = 0.0e20_wp
! TEST TWO  
!  blx(1:n) = -1.0_wp
!  bux(1:n) = +1.0e20_wp
! TEST THREE  
!  blx(1:n) = -1.0_wp
!  bux(1:n) = +1.0_wp
blx(1) = 0.59_wp!44
bux(1) = 1.445_wp
blx(5) = -3.9e-1_wp
bux(5) = -0.200_wp
! ------------------
bux(6) = 10.0_wp !! working good
! bux(6) = 11.0_wp !! fast !! 45 iterations 
! optimum at x =    0.6889229829244405   1.8735055049705871   2.0682924682678787
! 4.6402277568643573  -0.2444283933594384   1.8736741088671192
! 2.1732349926983754E-06
! bux(6) = 12.0_wp !! stagnation !! 43882 iterations
! optimum at x =    0.5900000000000000   2.0808044760419366   2.3127303035245053
! 5.0477989247787995  -0.3900000000000000   6.3782210775222410
! 6.4906367799130378E-06

! TEST four fail  
!   blx(1:n) = -1.0_wp
!   bux(1:n) = -1.1_wp

!   blx(1) = 0.087_wp
!   blx(2) = 0.955_wp
!   blx(3) = 0.85_wp
!   blx(4) = 2.96_wp
!   blx(5) = 1.59_wp
!   blx(6) = 4.99_wp
!   iusrbox = 1



  Call nlls_setup_bounds(params, n, blx, bux, options, inform)
  if (inform%status/=0) then
    Write(*,*) 'ERROR: nlls_setup_bounds failed, status=', inform%status
    stop
  End if

  Write(*,*) 'Box description (x0 is not proj)'
  Write(*,*) 'iusrbox = ', params%iusrbox
  if (params%iusrbox/=0) Then
    Do i = 1, n
      Write(*,Fmt=99999) blx(i), x(i), bux(i)
    End Do
  else
    Write(*,*) 'No bounds or all bound where -/+infinity'
  End If
99999  Format (5X,3(Es13.6e2,2X))
  options%print_level = 1
  options%exact_second_derivatives = .false. !.true.
!  options%model = 1 ! GN
!  options%model = 2 ! (Quasi-)Newton
  options%model = 3 ! Hybrid
!  options%model = 4 ! Newton-tensor
  options%type_of_method = 1 ! TR
!  options%type_of_method = 2 ! Regularization
!  options%nlls_method = 1 ! Powell's dogleg method
!  options%nlls_method = 2 ! Adachi-Iwata-Nakatsukasa-Takeda method
  options%nlls_method = 3 ! Moré-Sorensen method
!  options%nlls_method = 4 ! Galahan (DTRS:type_of_method=1, DRQS:type_of_method=2)
!  options%use_ews_subproblem = .true.
!  options%regularization_term = 1.0e-2
!  options%regularization_power = 2.0
!  options%reg_order = -1.0
!  options%inner_method = 1 ! passed in as a base reg term 
  options%inner_method = 2 ! expanded NLLS is solved
!  options%inner_method = 3 ! implicit recursive call
  options%maxit = 2000000
  options%box_linesearch_type = 1
!   options%box_tr_test_step = .False. !.True.
!   options%box_wolfe_test_step = .True.
!   options%box_max_ntrfail = 10

  call cpu_time(tic)
  call nlls_solve(n,m,x,eval_r, eval_J, eval_HF, params, options, inform)
  if(inform%status.ne.0) then
     print *, "ral_nlls() returned with error flag ", inform%status
     stop
  endif
  call cpu_time(toc)
  Write(*,*) 'Solution: '
  if (params%iusrbox/=0) Then
    Do i = 1, n
      Write(*,Fmt=99999) blx(i), x(i), bux(i)
    End Do
  else
    Do i = 1, n
      Write(*,Fmt=99999) -1.0e-20_wp, x(i), 1.0e20_wp
    End Do
  End If
  ! Print result
  print *, "Objective Value at solution    = ", inform%obj
  print *, "Objective Gradient at solution = ", inform%norm_g
  print *, "Found a local optimum at x = ", x
  print *, "Took ", inform%iter, " iterations (LS:",inform%ls_step_iter,' PG:',inform%pg_step_iter,')'
  print *, "     ", inform%f_eval, " function evaluations (LS:",inform%f_eval_ls,' PG:',inform%f_eval_pg,')'
  print *, "     ", inform%g_eval, " gradient evaluations (LS:",inform%g_eval_ls,' PG:',inform%g_eval_pg,')'
  print *, "     ", inform%h_eval, " hessian evaluations"
! print *, "     ", toc-tic, " seconds"
end program lanczos_box
