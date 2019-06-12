module example_module

  use :: ral_nlls_double
  use :: ral_nlls_internal
  use :: ral_nlls_workspaces
  implicit none 

  type, extends( params_base_type ) :: user_type
     real(wp), allocatable :: x_values(:)
     real(wp), allocatable :: y_values(:)
     integer :: m
  end type user_type

contains
  
  
SUBROUTINE eval_F( status, n_dummy, m, X, f, params)

!  -------------------------------------------------------------------
!  eval_F, a subroutine for evaluating the function f at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: f
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: X
       class( params_base_type ), intent(inout) :: params
! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )
       Real (Kind=wp) :: ex
       integer :: i
! then, let's work this into the format we need
! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          do i = 1,m
            ! Avoid overflow
            ex = max(-70.0_wp, min(70.0_wp, X(1) * params%x_values(i) + X(2)))
            f(i) = params%y_values(i) - exp( ex )
          end do
       end select

       status = 0
       
!!$! let's use Powell's function for now....
!!$       f(1) = X(1) + 10.0 * X(2)
!!$       f(2) = sqrt(5.0) * (X(3) - X(4))
!!$       f(3) = ( X(2) - 2.0 * X(3) )**2
!!$       f(4) = sqrt(10.0) * ( X(1) - X(4) )**2
       
! end of subroutine eval_F
       
     END SUBROUTINE eval_F

     subroutine eval_F_error( status, n_dummy, m_dummy, X_dummy, f_dummy, params_dummy)
       ! a fake eval_f to flag an error 
       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: f_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: X_dummy
       class( params_base_type ), intent(inout) :: params_dummy

       status = -1
       
     end subroutine eval_F_error

     SUBROUTINE eval_J( status, n_dummy, m, X, J, params)

!  -------------------------------------------------------------------
!  eval_J, a subroutine for evaluating the Jacobian J at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: J
       REAL ( wp ), DIMENSION( * ),INTENT( IN ) :: X
       class( params_base_type ), intent(inout) :: params

! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i

       ! let's work this into the format we need
       ! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          do i = 1,m
             J(i) =  - params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) )
             J(m + i) = - exp( X(1) * params%x_values(i) + X(2) )
          end do
       end select
       
       status = 0

! end of subroutine eval_J

!!$       ! initialize to zeros...
!!$       J(1:4,1:4) = 0.0
!!$       
!!$       ! enter non-zeros values
!!$       J(1,1) = 1.0
!!$       J(1,2) = 10.0
!!$       J(2,3) = sqrt(5.0)
!!$       J(2,4) = -sqrt(5.0)
!!$       J(3,2) = 2.0 * (X(2) - 2.0 * X(3))
!!$       J(3,3) = -4.0 * (X(2) - 2.0 * X(3)) 
!!$       J(4,1) = sqrt(10.0) * 2.0 * (X(1) - X(4))
!!$       J(4,4) = - sqrt(10.0) * 2.0 * (X(1) - X(4))

     END SUBROUTINE eval_J

     SUBROUTINE eval_J_c( status, n_dummy, m, X, J, params)

!  -------------------------------------------------------------------
!  eval_J, a subroutine for evaluating the Jacobian J at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: J
       REAL ( wp ), DIMENSION( * ),INTENT( IN ) :: X
       class( params_base_type ), intent(inout) :: params

! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i, index

       ! let's work this into the format we need
       ! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          index = 0
          ! use c-based formatting
          do i = 1, m
             J(index + 1) = - params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) )
             J(index + 2) = - exp( X(1) * params%x_values(i) + X(2) )
             index = index + 2
          end do
       end select
       
       status = 0

! end of subroutine eval_J

!!$       ! initialize to zeros...
!!$       J(1:4,1:4) = 0.0
!!$       
!!$       ! enter non-zeros values
!!$       J(1,1) = 1.0
!!$       J(1,2) = 10.0
!!$       J(2,3) = sqrt(5.0)
!!$       J(2,4) = -sqrt(5.0)
!!$       J(3,2) = 2.0 * (X(2) - 2.0 * X(3))
!!$       J(3,3) = -4.0 * (X(2) - 2.0 * X(3)) 
!!$       J(4,1) = sqrt(10.0) * 2.0 * (X(1) - X(4))
!!$       J(4,4) = - sqrt(10.0) * 2.0 * (X(1) - X(4))

     END SUBROUTINE eval_J_c


     
     subroutine eval_J_error( status, n_dummy, m_dummy, X_dummy, J_dummy, params_dummy)
       ! a fake eval_J to flag an error 
       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: J_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: X_dummy
       class( params_base_type ), intent(inout) :: params_dummy

       status = -1
       
     end subroutine eval_J_error

     SUBROUTINE eval_H( status, n_dummy, m, X, f, h, params)

!  -------------------------------------------------------------------
!  eval_H, a subroutine for evaluating the second derivative hessian terms
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: f
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: h
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: X
       class( params_base_type ), intent(inout) :: params
! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i

! then, let's work this into the format we need
! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          ! evaluate 
          ! HF = \sum_{i=1}^m F_i H_i
          h(1:4) = 0.0
          do i = 1, m
             h(1) = &
                  h(1) + f(i)* ( & 
                  - (params%x_values(i)**2) * exp( X(1) * params%x_values(i) + X(2) ) &
                  )
             h(2) = &
                  h(2) + f(i)* ( &
                  - params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) ) &
                  )
             h(4) = &
                  h(4) + f(i)* ( &
                  -  exp( X(1) * params%x_values(i) + X(2) ) &
                  )
          end do
          h(3) = h(2)
       end select

       status = 0
       
!!$! let's use Powell's function for now....
!!$       f(1) = X(1) + 10.0 * X(2)
!!$       f(2) = sqrt(5.0) * (X(3) - X(4))
!!$       f(3) = ( X(2) - 2.0 * X(3) )**2
!!$       f(4) = sqrt(10.0) * ( X(1) - X(4) )**2
       
! end of subroutine eval_F
       
     END SUBROUTINE eval_H

     subroutine eval_H_error( status, n_dummy, m_dummy, X_dummy, f_dummy, h_dummy, params_dummy)

!  -------------------------------------------------------------------
!  a fake eval_H for flagging an error
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n_dummy, m_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: f_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: h_dummy
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: X_dummy
       class( params_base_type ), intent(inout) :: params_dummy

       status = -1
       
     end subroutine eval_H_error

     subroutine eval_HP ( status, n, m, x, y, hp, params )
       
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: status
       INTEGER, INTENT( IN ) :: n, m 
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: x
       REAL ( wp ), DIMENSION( * ),INTENT( IN )  :: y
       REAL ( wp ), DIMENSION( * ),INTENT( OUT ) :: hp
       class( params_base_type ), intent(inout) :: params

       ! does nothing for now...

       integer :: i

       ! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)

          hp(1:n*m) = 0.0
          do i = 1, m ! loop over the columns
             ! need to put H(x)*y in each row
             hp( n*(i-1) + 1 ) = &
                  y(1)* (- (params%x_values(i)**2) * exp( X(1) * params%x_values(i) + X(2) ) ) + &
                  y(2)* (- params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) ) )
             hp( n*(i-1) + 2 ) = &
                  y(1)* (- params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) ) ) + &
                  y(2)* (-  exp( X(1) * params%x_values(i) + X(2) ) )
          end do
 
       end select
       
       status = 0
       
     end subroutine eval_HP

     subroutine generate_data_example(params)
       
       type ( user_type ), intent(out) :: params
       
       params%m = 67
       allocate(params%x_values(params%m))
       allocate(params%y_values(params%m))
              
       ! First, let's get the data
       ! Generated with the code
       !  randn('seed', 23497);
       !   m = 0.3;
       !   c = 0.1;
       !   x_data = [0:0.075:5];
       !   y = exp(m * x_data + c);
       !   noise = randn(size(x_data)) * 0.2;
       !   y_data = y + noise;
       ! (c.f. https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/curve_fitting.cc)
       params%x_values = (/ 0.0, &
   0.075000000000000, &
   0.150000000000000, &
   0.225000000000000, &
   0.300000000000000, &
   0.375000000000000, &
   0.450000000000000, &
   0.525000000000000, &
   0.600000000000000, &
   0.675000000000000, &
   0.750000000000000, &
   0.825000000000000, &
   0.900000000000000, &
   0.975000000000000, &
   1.050000000000000, &
   1.125000000000000, &
   1.200000000000000, &
   1.275000000000000, &
   1.350000000000000, &
   1.425000000000000, &
   1.500000000000000, &
   1.575000000000000, &
   1.650000000000000, &
   1.725000000000000, &
   1.800000000000000, &
   1.875000000000000, &
   1.950000000000000, &
   2.025000000000000, &
   2.100000000000000, &
   2.175000000000000, &
   2.250000000000000, &
   2.325000000000000, &
   2.400000000000000, &
   2.475000000000000, &
   2.550000000000000, &
   2.625000000000000, &
   2.700000000000000, &
   2.775000000000000, &
   2.850000000000000, &
   2.925000000000000, &
   3.000000000000000, &
   3.075000000000000, &
   3.150000000000000, &
   3.225000000000001, &
   3.300000000000000, &
   3.375000000000000, &
   3.450000000000000, &
   3.525000000000000, &
   3.600000000000001, &
   3.675000000000000, &
   3.750000000000000, &
   3.825000000000000, &
   3.900000000000000, &
   3.975000000000000, &
   4.050000000000001, &
   4.125000000000000, &
   4.200000000000000, &
   4.275000000000000, &
   4.350000000000001, &
   4.425000000000000, &
   4.500000000000000, &
   4.575000000000000, &
   4.650000000000000, &
   4.725000000000001, &
   4.800000000000000, &
   4.875000000000000, &
   4.950000000000000 /)
       
       params%y_values = (/ 0.907946872110432, &
   1.199579396036134, &
   1.060092431384317, &
   1.298370500472354, &
   0.952768858414788, &
   1.209665290655204, &
   1.256912538155493, &
   1.163922146095987, &
   1.004877938808100, &
   1.205944250961060, &
   0.952693297695969, &
   1.449662692280761, &
   1.402015259144406, &
   1.378094012325746, &
   1.560882147577552, &
   1.437185539058121, &
   1.559853079888265, &
   1.877814947316832, &
   1.818781749024682, &
   1.375546045112591, &
   1.233967904388409, &
   1.887793124397751, &
   1.610237096463521, &
   1.787032484792262, &
   1.850015127982676, &
   2.120553361509177, &
   1.942913663511919, &
   2.106517132599766, &
   2.271787117356578, &
   1.727554346001754, &
   2.002909500898113, &
   1.975837413903495, &
   2.337446525801909, &
   1.960190841677278, &
   2.447097025572309, &
   2.161663720225506, &
   2.748798529374621, &
   2.507814238594416, &
   2.423769408403069, &
   2.578119353028746, &
   2.460310096221557, &
   2.638362783992324, &
   2.765540456237868, &
   2.837165966564409, &
   3.179711963042789, &
   3.245315453091675, &
   3.289631922410174, &
   3.360995198615834, &
   3.470489725998371, &
   3.169513520153466, &
   3.363740517933189, &
   3.665288099084969, &
   3.620334359722351, &
   4.018911445550667, &
   3.512715166706162, &
   3.874661411575566, &
   4.197746303653517, &
   3.703511523106007, &
   4.076351488309604, &
   4.056340365649961, &
   4.297751562451419, &
   4.373076571153739, &
   4.577093065941748, &
   4.856619059058190, &
   4.927350280596274, &
   4.703122139742729, &
   4.870205182453842 /)
       
     end subroutine generate_data_example
     
     subroutine reset_default_options(options)
       type( nlls_options ), intent(inout) :: options
       
       type( nlls_options ) :: default_options

!       options%error = default_options%error
!       options%out = default_options%out
!       options%print_level = default_options%print_level
       options%print_options = default_options%print_options
       options%print_header = default_options%print_header 
       options%maxit = default_options%maxit
       options%model = default_options%model
       options%type_of_method = default_options%type_of_method
       options%nlls_method = default_options%nlls_method
       options%lls_solver = default_options%lls_solver
       options%stop_g_absolute = default_options%stop_g_absolute
       options%stop_g_relative = default_options%stop_g_relative
       options%stop_f_absolute = default_options%stop_f_absolute
       options%stop_f_relative = default_options%stop_f_relative
       options%stop_s = default_options%stop_s       
       options%relative_tr_radius = default_options%relative_tr_radius
       options%initial_radius_scale = default_options%initial_radius_scale
       options%initial_radius = default_options%initial_radius
       options%base_regularization = default_options%base_regularization
       options%regularization = default_options%regularization
       options%regularization_term = default_options%regularization_term
       options%regularization_power = default_options%regularization_power
       options%maximum_radius = default_options%maximum_radius
       options%eta_successful = default_options%eta_successful
       options%eta_success_but_reduce = default_options%eta_success_but_reduce
       options%eta_very_successful = default_options%eta_very_successful
       options%eta_too_successful = default_options%eta_too_successful
       options%radius_increase = default_options%radius_increase
       options%radius_reduce = default_options%radius_reduce
       options%radius_reduce_max = default_options%radius_reduce_max
       options%tr_update_strategy = default_options%tr_update_strategy
       options%hybrid_switch = default_options%hybrid_switch
       options%exact_second_derivatives = default_options%exact_second_derivatives
       options%subproblem_eig_fact = default_options%subproblem_eig_fact
       options%scale = default_options%scale
       options%scale_max = default_options%scale_max
       options%scale_min = default_options%scale_min
       options%scale_trim_max = default_options%scale_trim_max
       options%scale_trim_min = default_options%scale_trim_min
       options%scale_require_increase = default_options%scale_require_increase
       options%setup_workspaces = default_options%setup_workspaces
       options%remove_workspaces = default_options%remove_workspaces
       options%more_sorensen_maxits = default_options%more_sorensen_maxits
       options%more_sorensen_shift = default_options%more_sorensen_shift
       options%more_sorensen_tiny = default_options%more_sorensen_tiny
       options%more_sorensen_tol = default_options%more_sorensen_tol
       options%hybrid_tol = default_options%hybrid_tol
       options%hybrid_switch_its = default_options%hybrid_switch_its
       options%reg_order = default_options%reg_order
       options%inner_method = default_options%inner_method
       options%output_progress_vectors = default_options%output_progress_vectors
       options%update_lower_order = default_options%update_lower_order
       
     end subroutine reset_default_options

     subroutine solve_basic(X,params,options,inform)
      
       real(wp), intent(out) :: X(:)
       type( user_type ), intent(inout) :: params
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform

       integer :: n, m
       
       n = 2 
       m = 67
       
       X(1) = 1.0
       X(2) = 2.0
       

       call nlls_solve(n, m, X,                         &
                   eval_F, eval_J, eval_H, params,  &
                   options, inform )       
       
     end subroutine solve_basic

     subroutine dogleg_tests(options,fails)
       
       type( nlls_options ), intent(inout) :: options       
       integer, intent(out) :: fails
       
       real(wp), allocatable :: J(:), hf(:), f(:), g(:), d(:)
       real(wp) :: Delta, normd
       type( nlls_inform ) :: inform
       type( nlls_workspace ) :: w
       type( nlls_workspace ), Target :: iw

       integer :: n,m

       fails = 0
!      Link the inner_workspace to the main workspace
       w%iw_ptr => iw
!      Self reference for inner workspace so recursive call does not fail
       iw%iw_ptr => iw

       options%scale = 0 
       options%print_level = 3
       
       !! dogleg 
       options%nlls_method = 1
       options%model = 5
       n = 2
       m = 3
       allocate(J(m*n), hf(n*n), f(m), g(n), d(n))
       call setup_workspaces(w,n,m,options,inform) 
       Delta = 10.0_wp
     
       ! first, hit the 'method not supported' error
       options%model = 27
       J = 1.0_wp
       hf = 0.0_wp
       f = 1.0_wp
       g = 1.0_wp
       call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%calculate_step_ws%dogleg_ws)
       if (inform%status .ne. NLLS_ERROR_DOGLEG_MODEL) then
          write(*,*) 'Error: unsupported model allowed in dogleg'
          fails = fails + 1
       end if
       inform%status = 0
       
       options%model = 1
       J  = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
       f = 1.0_wp
       hf = 0.0_wp
       g  = 1.0_wp
       ! now, get ||d_gn|| <= Delta
       Delta = 6.0_wp
       call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%calculate_step_ws%dogleg_ws)
       if (inform%status .ne. 0) then
          write(*,*) 'Error: unexpected error in dogleg'
          fails = fails + 1
       end if

       ! now set delta so that || alpha * d_sd || >= Delta
       Delta = 0.5_wp
       call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%calculate_step_ws%dogleg_ws)
       if (inform%status .ne. 0) then
          write(*,*) 'Error: unexpected error in dogleg'
          fails = fails + 1
       end if

       ! now get the guys in the middle...
       Delta = 2.5_wp
       call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%calculate_step_ws%dogleg_ws)
       if (inform%status .ne. 0) then
          write(*,*) 'Error: unexpected error in dogleg'
          fails = fails + 1
          inform%status = 0
       end if
     
       call nlls_finalize(w,options)
     
       call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%calculate_step_ws%dogleg_ws)
       if (inform%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
          write(*,*) 'Error: workspace error not flagged when workspaces not setup'
          fails = fails + 1
       end if

       call reset_default_options(options)
       
     end subroutine dogleg_tests

     subroutine generate_scaling_tests(options, fails)
       type(nlls_options),intent(inout) :: options 
       integer, intent(out) :: fails

       real(wp), allocatable :: J(:), A(:,:), scale_extra(:), scale(:)
       integer :: n,m
       type( nlls_workspace) :: w
       type( nlls_workspace), Target :: iw
       type( nlls_inform ) :: inform
       
       fails = 0
       w%iw_ptr => iw
       iw%iw_ptr => iw

       options%scale = 4
       options%nlls_method = 3
       n = 2
       m = 3 
       allocate(J(m*n),A(n,n),scale_extra(n),scale(n))
       call setup_workspaces(w,n,m,options,inform)
     
       J = 0.0_wp
       J(1) = 1e15
       J(4) = 1e-15
       A(1,1) = 1.0_wp
       A(2,1) = 0.0_wp
       A(1,2) = 0.0_wp
       A(2,2) = 1.0_wp

       scale = 1.0_wp
       scale_extra = 0.0_wp 

       !** scale = 1 **
       options%scale= 1     
       call generate_scaling(J,A,n,m,scale,scale_extra,& 
            w%calculate_step_ws%generate_scaling_ws, &
            options,inform)
       if (inform%status .ne. 0 ) then
          write(*,*) 'Error: unexpected error in generate_scaling when scale = 1'
          write(*,*) 'status = ', inform%status,' returned.'
          fails = fails + 1
          inform%status = 0 
       end if

       !** scale = 2 **
       options%scale = 2
       call generate_scaling(J,A,n,m,scale,scale_extra,& 
            w%calculate_step_ws%generate_scaling_ws, &
            options,inform)
       if (inform%status .ne. 0 ) then
          write(*,*) 'Error: unexpected error in generate_scaling when scale = 2'
          write(*,*) 'status = ', inform%status,' returned.'
          fails = fails + 1
          inform%status = 0 
       end if

       !** scale = 3 **
       options%scale = 3
       call generate_scaling(J,A,n,m,scale,scale_extra,& 
            w%calculate_step_ws%generate_scaling_ws, &
            options,inform)
       if (inform%status .ne. 0 ) then
          write(*,*) 'Error: unexpected error in generate_scaling when scale = 2'
          write(*,*) 'status = ', inform%status,' returned.'
          fails = fails + 1
          inform%status = 0 
       end if

     !** scale undefined
     options%scale = 786
     call generate_scaling(J,A,n,m,scale,scale_extra,& 
          w%calculate_step_ws%generate_scaling_ws, &
          options,inform)
     if (inform%status .ne. NLLS_ERROR_BAD_SCALING ) then
        write(*,*) 'Error: expected error in generate_scaling when passing undefined scaling'
        write(*,*) 'status = ', inform%status,' returned.'
        fails = fails + 1
        inform%status = 0 
     end if
     inform%status = 0
     
     ! now, let's test the non-default modes
     ! first, set scale_require_increase to T
     options%scale = 1
     options%scale_require_increase = .true.
     call generate_scaling(J,A,n,m,scale,scale_extra,& 
          w%calculate_step_ws%generate_scaling_ws, &
          options,inform)
     if (inform%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', inform%status,' returned.'
        fails = fails + 1
        inform%status = 0 
     end if
     options%scale_require_increase = .false.

     ! first, set scale_trim_min to T
     options%scale_trim_min = .true.
     call generate_scaling(J,A,n,m,scale,scale_extra,& 
          w%calculate_step_ws%generate_scaling_ws, &
          options,inform)
     if (inform%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', inform%status,' returned.'
        fails = fails + 1
        inform%status = 0 
     end if
     options%scale_trim_min = .false.

     ! first, set scale_trim_max to T
     options%scale_trim_max = .false.
     call generate_scaling(J,A,n,m,scale,scale_extra,& 
          w%calculate_step_ws%generate_scaling_ws, &
          options,inform)
     if (inform%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', inform%status,' returned.'
        fails = fails + 1
        inform%status = 0 
     end if
     options%scale_trim_max = .true.


     call nlls_finalize(w,options)
     call generate_scaling(J,A,n,m,scale,scale_extra,& 
          w%calculate_step_ws%generate_scaling_ws, &
          options,inform)
     if (inform%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        write(*,*) '(generate_scaling)'
        fails = fails + 1
     end if

     call reset_default_options(options)

   end subroutine generate_scaling_tests

   subroutine aint_tr_tests(options,fails)
     
     type( nlls_options), intent(inout) :: options
     integer, intent(out) :: fails
     
     REAL(wp), allocatable :: J(:), A(:,:), hf(:), f(:), v(:), X(:), d(:)
     real(wp) :: Delta, normd
     integer  :: n, m
     TYPE( nlls_inform ) :: inform
     type( nlls_workspace ) :: w
     type( nlls_workspace), Target :: iw
      
     fails = 0
     w%iw_ptr => iw
     iw%iw_ptr => iw
     options%nlls_method = 2
     
     n = 2
     m = 3

     allocate(J(n*m), A(n,n), f(m), v(n), X(n), hf(n*n),d(n))
     
     J = [ 1.0, 1.0, 2.0, 2.0, 3.0, 4.0 ]
     A = reshape([6.0, 13.0, 13.0, 29.0],[2, 2])
     f = [2.0, 3.0, 4.0]
     v = [13.0, 29.0]
     hf = 0.0_wp
     Delta = 1.0_wp
          
     call setup_workspaces(w,n,m,options,inform) 
     call aint_tr(J,A,f,X,v,hf,n,m,Delta,d,normd,options,inform,& 
          w%calculate_step_ws%aint_tr_ws)
     if (inform%status .ne. 0) then 
        write(*,*) 'Error: aint_tr test failed'
        fails = fails + 1
     end if
     
     call nlls_finalize(w,options)
     
     call aint_tr(J,A,f,X,v,hf,n,m,Delta,d,normd,options,inform,& 
          w%calculate_step_ws%aint_tr_ws)
     if (inform%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        write(*,*) '(aint_tr)'
        fails = fails + 1
     end if
     call reset_default_options(options)

   end subroutine aint_tr_tests

   subroutine more_sorensen_tests(options,fails)
     type( nlls_options), intent(inout) :: options
     integer, intent(out) :: fails
     
     
     real(wp), allocatable :: g(:), d(:), A(:,:)
     real(wp) :: Delta, normd
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     integer :: n,m
     type( nlls_workspace ), Target :: iw
      
     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw
     
     options%nlls_method = 3
     
     n = 2
     m = 3
     allocate(A(n,n), g(n), d(n))
     call setup_workspaces(work,n,m,options,status) 
     Delta = 10.0_wp
     
     ! regular case...
     A(1,1) = 0.2_wp
     A(2,1) = 0.3_wp
     A(1,2) = A(2,1)
     A(2,2) = 0.4_wp
     g = 1.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen'
        fails = fails + 1
        status%status = 0
     end if

     ! non spd matrix, with failure
     options%more_sorensen_shift = -1000.0_wp
     g = 1.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. NLLS_ERROR_MS_TOO_MANY_SHIFTS) then
        write(*,*) 'Error: MS too many shifts test passed, when fail expected'
        fails = fails + 1
     end if
     status%status = 0
     options%more_sorensen_shift = 1e-13

     ! look for nd /=  Delta with a non-zero shift?
     g = 1.0_wp
     Delta =  10.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with non-zero shift'
        write(*,*) 'status = ', status%status, ' returned'
        fails = fails + 1
        status%status = 0
     end if

     ! now look for nd =  Delta with a non-zero shift?
     g = 1.0_wp
     options%more_sorensen_tiny = 0.01_wp
     Delta =  0.2055_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with non-zero shift'
        write(*,*) 'status = ', status%status, ' returned'
        fails = fails + 1
        status%status = 0
     end if
     call reset_default_options(options)
     options%nlls_method = 3
!     options%more_sorensen_tiny = beta
     ! *todo*

     ! now take nd > Delta
     d = 1.0_wp
     Delta = 3.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with nd > Delta'
        fails = fails + 1
        status%status = 0
     end if


     ! get to max_its...
     options%more_sorensen_maxits = 1     
     g = 1.0_wp
     Delta = 3.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. NLLS_ERROR_MS_MAXITS) then
        write(*,*) 'Error: Expected maximum iterations error in more_sorensen'
        fails = fails + 1
     end if
     status%status = 0
     options%more_sorensen_maxits = 10
     
     call nlls_finalize(work,options)
     call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1
     end if

     call reset_default_options(options)
     
   end subroutine more_sorensen_tests

   subroutine trust_region_subproblem_tests(options,fails)

     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: g(:), d(:), A(:,:)
     real(wp) :: Delta, normd
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     integer :: n,m
     integer :: problem, method
     character (len=40) :: problem_name
     character (len=40) :: method_name
     integer :: num_successful_steps
     
     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw

     options%out = 6
     options%print_level = 0
     
     n = 4
     m = 5
     allocate(A(n,n), g(n), d(n))

     do problem = 1, 6
        select case (problem)
        case (1,2) ! Convex, within the tr radius
           g = 1.0_wp
           A = 0.0_wp
           A(1,1) = 1.0_wp
           A(2,2) = 2.0_wp
           A(3,3) = 3.0_wp
           A(4,4) = 4.0_wp

           if (problem == 1) then
              Delta = 1.0_wp ! point lies on the tr radius
           elseif (problem == 2) then 
              Delta = 3.0_wp ! point lies in the tr radius
           end if
           
           write(problem_name,'(A,ES12.4)') '7.3.1.1, Delta = ',Delta
        case (3,4)  ! Nonconvex
           g = 1.0_wp
           A = 0.0_wp           
           A(1,1) = -2.0_wp
           A(2,2) = -1.0_wp
           A(3,3) =  0.0_wp
           A(4,4) =  1.0_wp
           if (problem == 3) then
              Delta = 1.0_wp ! point lies on the tr radius
           elseif (problem == 4) then 
              Delta = 3.0_wp ! point lies in the tr radius
           end if
           
           write(problem_name,'(A,ES12.4)') '7.3.1.2, Delta = ',Delta
        case (5,6) ! hard case
           g = 1.0_wp
           g(1) = 0.0_wp
           A = 0.0_wp           
           A(1,1) = -2.0_wp
           A(2,2) = -1.0_wp
           A(3,3) =  0.0_wp
           A(4,4) =  1.0_wp
           if (problem == 5) then
              Delta = 1.0_wp ! point lies on the tr radius
           elseif (problem == 6) then 
              Delta = 1.5_wp ! point lies in the tr radius
           end if
           write(problem_name,'(A,ES12.4)') '7.3.1.3, Delta = ',Delta
        end select
        do method = 1,3 ! now, let's loop through the methods available...
           select case (method)
           case (1) ! more sorensen, eigenvalues
              method_name = 'more sorensen, eigenvalues'
              options%nlls_method = 3 
              options%use_ews_subproblem = .true.
              call setup_workspaces(work,n,m,options,status)
              call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
                   work%calculate_step_ws%more_sorensen_ws)
           case (2) ! dltr
              method_name = 'dltr'
              options%nlls_method = 4
              num_successful_steps = 0 
              call setup_workspaces(work,n,m,options,status)
              call solve_galahad(A,g,n,m,Delta,num_successful_steps,& 
                   d,normd,2.0_wp,options,status,&
                   work%calculate_step_ws%solve_galahad_ws )
           case (3) ! more sorensen, no eigenvalues
              method_name = 'more_sorensen, no eigenvalues'
              options%nlls_method = 3
              options%use_ews_subproblem = .false.
              call setup_workspaces(work,n,m,options,status)
              call more_sorensen(A,g,n,m,Delta,d,normd,options,status,& 
                   work%calculate_step_ws%more_sorensen_ws)              
           end select
           if (options%print_level > 0 ) then 
              write(*,*) 'method = ', method_name, ' problem = ', problem_name, 'normd = ', normd
           end if
           if ( (status%status .ne. 0)) then
              write(*,*) 'Error: unexpected error in ', method_name
              write(*,*) '(status = ',status%status,')'
              write(*,*) 'TR Book Example ',problem_name
              fails = fails + 1
              status%status = 0
           elseif ( normd - Delta > 1e-3 ) then
              write(*,*) 'Error: answer returned outside the TR Radius'
              write(*,*) 'TR Book Example ', trim(problem_name),' using method ',method_name
              write(*,*) 'Delta = ', Delta, '||d|| = ', normd
              fails = fails + 1
              status%status = 0
           end if
        end do       
     end do

     options%out = 17
     
   end subroutine trust_region_subproblem_tests
   
   subroutine evaluate_model_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: f(:), J(:), hf(:), X(:), Xnew(:), d(:)
     real(wp) :: md, md_gn
     integer :: m, n
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     
     real(wp) :: one = 1.0_wp

     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw
     n = 2
     m = 4
     
     allocate(f(m), J(n*m), hf(n*n), X(n), Xnew(n), d(n))
     f = one
     J = one
     hf = one 
     X = one
     Xnew = one
          
     call evaluate_model(f,J,hf,X,Xnew,d,md,md_gn,m,n,options,status,& 
          work%calculate_step_ws%evaluate_model_ws)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1
     end if
     
     
     call reset_default_options(options)
     
   end subroutine evaluate_model_tests

   subroutine solve_galahad_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: g(:), d(:), A(:,:)
     real(wp) :: Delta, normd
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     integer :: n,m, i, num_successful_steps
     character (len=5) :: testname
     
     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw

     options%nlls_method = 4
     
     n = 2
     m = 5
     
     num_successful_steps = 0 
     
     allocate(A(n,n), g(n), d(n))

     do i = 1,2
        options%type_of_method = i
        if (i == 1) then 
           testname = 'DTRS '
        elseif (i == 2) then 
           testname = 'DRQS '
        end if
        
        call setup_workspaces(work,n,m,options,status) 
        
        A(1,1) = 10.0
        A(1,2) = 2.0
        A(2,1) = 2.0
        A(2,2) = 10.0
        g = [-7.4, -28.9]
        Delta = 0.02_wp
        call solve_galahad(A,g,n,m,Delta,num_successful_steps,& 
             d,normd,2.0_wp,options,status,&
             work%calculate_step_ws%solve_galahad_ws )
        if ( status%status .ne. 0 ) then
           write(*,*) testname,'test failed, status = ', status%status
           fails = fails + 1
        end if
        
        if (i == 1) then
           ! check result lies within the trust region
           if ( abs(dot_product(d,d) - Delta**2) > 1e-12 ) then
              write(*,*) testname,'failed'
              write(*,*) 'Delta = ', Delta, '||d|| = ', dot_product(d,d)
              fails = fails + 1
           end if
        end if

        
        Delta = -100.0_wp
        call solve_galahad(A,g,n,m,Delta,num_successful_steps,& 
             d,normd,2.0_wp,options,status,&
             work%calculate_step_ws%solve_galahad_ws )
        if ( status%status .ne. NLLS_ERROR_FROM_EXTERNAL ) then
           write(*,*) testname,'test failed, expected status = ', NLLS_ERROR_FROM_EXTERNAL
           write(*,*) ' but got status = ', status%status
           fails = fails + 1
        end if
        status%status = 0

        call nlls_finalize(work,options)
                
     end do

     call reset_default_options(options)
     
   end subroutine solve_galahad_tests

   subroutine solve_newton_tensor_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: J(:), f(:), X(:), d(:)
     real(wp) :: Delta, md
     integer :: n, m, num_successful_steps
     type( params_base_type ) :: params
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( tenJ_type ), Target :: tenJ
     type( tenJ_type ), Pointer :: tenJ_pointer
     type( NLLS_workspace ), Target :: inner_workspace
     real(wp) :: one = 1.0_wp
     
     fails = 0
     tenJ_pointer => tenJ
     work%iw_ptr => inner_workspace
     inner_workspace%iw_ptr => inner_workspace
     
     n = 3 
     m = 5
     allocate(J(n*m),f(m),X(n),d(n))
     J = one
     f = one
     X = one

     num_successful_steps = 0
     
     call solve_newton_tensor(J, f, eval_H, X, n, m, Delta, num_successful_steps, & 
                                    d, md, params, options, status, & 
                                    work%calculate_step_ws%solve_newton_tensor_ws,&
                                    tenJ_pointer, inner_workspace)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1
     end if
         
     
     call reset_default_options(options)
     
   end subroutine solve_newton_tensor_tests

   subroutine all_eig_symm_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), ew(:), ev(:,:)
     integer :: n
     type( nlls_workspace ) :: work
     type( nlls_inform ) :: status
     
     real(wp) :: one = 1.0_wp

     fails = 0

     n = 2
     allocate(A(n,n), ew(n), ev(n,n))
     A = one

     call all_eig_symm(A,n,ew,ev, & 
             work%calculate_step_ws%solve_galahad_ws%all_eig_symm_ws, status)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1
     end if
          
     call reset_default_options(options)
     
   end subroutine all_eig_symm_tests

   subroutine solve_LLS_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: J(:), f(:), d(:), Jd(:)
     real(wp) :: normerror
     integer :: n,m 
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     type( nlls_inform ) :: status

     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw

     options%nlls_method = 1 ! dogleg

     n = 2
     m = 5

     call setup_workspaces(work,n,m,options,status) 

     allocate(J(n*m), f(m), d(n), Jd(m))
     J = [ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, & 
           6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp ]
     f = [ 7.0_wp, 9.0_wp, 11.0_wp, 13.0_wp, 15.0_wp ]

     call solve_LLS(J,f,n,m,d,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if ( status%status .ne. 0 ) then 
        write(*,*) 'solve_LLS test failed: wrong error message returned'
        write(*,*) 'status = ', status%status, " (expected ",NLLS_ERROR_FROM_EXTERNAL,")"
        fails = fails + 1
     end if
     ! check answer
     call mult_J(J,n,m,d,Jd)
     normerror = norm2(Jd + f)
     if ( normerror > 1.0e-12_wp ) then
        ! wrong answer, as data chosen to fit
        write(*,*) 'solve_LLS test failed: wrong solution returned'
        write(*,*) '||Jd - f|| = ', normerror
        fails = fails + 1
     end if
     
     n = 100 
     m = 20
     deallocate(J,f,Jd,d)
     allocate(J(n*m), f(m), d(n)) ! f has wrong size
     
     call setup_workspaces(work,n,m,options,status) 

     J = 1.0_wp
     f = 1.0_wp
     call solve_LLS(J,f,n,m,d,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if ( status%status .ne. NLLS_ERROR_FROM_EXTERNAL ) then 
        write(*,*) 'solve_LLS test failed: wrong error message returned'
        write(*,*) 'status = ', status%status, " (expected ",NLLS_ERROR_FROM_EXTERNAL,")"
        fails = fails + 1
     end if
     status%status = 0

     call nlls_finalize(work, options)
     
     call solve_LLS(J,f,n,m,d,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1
     end if

     
     call reset_default_options(options)
     
   end subroutine solve_LLS_tests


   subroutine findbeta_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: a(:), b(:)
     real(wp) :: Delta, beta
     type( nlls_inform ) :: status
     integer :: n
     
     fails = 0 

     n = 3
     allocate(a(n), b(n))
     a = [ 1.0, 2.0, 3.0 ]
     b = [ 2.0, 1.0, 1.0 ]

     call findbeta(a,b,10.0_wp,beta,status)

     if (status%status .ne. 0) then
        write(*,*) 'error -- findbeta did not work: info /= 0'
        fails = fails + 1
     else if ( ( norm2( a + beta * b ) - 10.0_wp ) > 1e-12 ) then
        write(*,*) 'error -- findbeta did not work'
        write(*,*) '|| x + beta y|| = ', norm2( (a + beta * b)-10.0_wp)
        fails = fails + 1
     end if
     
     deallocate(a,b)
     
     n = 2
     allocate(a(n),b(n))!,z(n))
     
     a(1) = 1e25_wp
     a(2) = 0.0_wp
     b(1) = 0.0_wp
     b(2) = 1.0_wp
     Delta = 1e-8_wp
     beta = 0.0_wp

     call findbeta(a,b,Delta,beta,status)

     if (status%status .ne. NLLS_ERROR_FIND_BETA) then
        write(*,*) 'Expected an error from findbeta: info =', status%status
        write(*,*) 'beta returned = ', beta
        write(*,*) '|| x + beta y|| = ', norm2( (a + beta * b) )
        fails = fails + 1
     end if

     call reset_default_options(options)
     
   end subroutine findbeta_tests

   subroutine calculate_rho_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp) :: normf, normfnew, md, rho

     fails = 0 

     normf = 2.0_wp
     normfnew = 1.0_wp
     md = 1.5_wp
     call calculate_rho(normf, normfnew, md, rho,options)
     if ( abs(rho - 3.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 3.0, got ', rho
        fails = fails + 1
     end if
     
     ! now, let's check one is returned if alpha = beta
     normfnew = 2.0_wp
     call calculate_rho(normf, normfnew, md, rho,options)
     if (abs(rho - 1.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 1.0, got ', rho
        fails = fails + 1
     end if
     normfnew = 1.0_wp

     ! finally, check that 1 is returned if denominator = 0
     md = 2.0_wp
     call calculate_rho(normf, normfnew, md, rho,options)
     if (abs(rho - 1.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 1.0, got ', rho
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine calculate_rho_tests

   subroutine update_trust_region_radius_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp) :: rho
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     
     fails = 0 
     work%iw_ptr => iw
     iw%iw_ptr => iw

     call setup_workspaces(work,2,2,options,status) 
     work%Delta = 100.0_wp ! Delta
     work%tr_nu = 2.0_wp ! nu
     work%tr_p = 3 ! p
     ! rho = rho

     options%tr_update_strategy = 1
     ! let's go through the options
     
     options%eta_success_but_reduce = 0.25_wp
     options%eta_very_successful = 0.75_wp
     options%eta_too_successful = 2.0_wp

     ! check if rho reduced...
     rho = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( work%Delta >= 100.0_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp
     
     ! check if rho stays the same...
     rho = (options%eta_success_but_reduce + options%eta_very_successful) / 2
     call update_trust_region_radius(rho,options,status,work)
     if ( abs(work%Delta - 100.0_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: Delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp

     ! check if rho increases...
     rho = (options%eta_very_successful + options%eta_too_successful) / 2
     work%norm_S_d = 100.0_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( work%Delta <= 100.0_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp

     
     ! check if rho stays the same because too successful...
     rho = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( abs(work%Delta - 100.0_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp

     ! now check for NaNs...HOW to do this in a non-compiler dependent way!?!?

     !! now, let's check the other option....
     options%tr_update_strategy = 2
     
     ! check if rho increases...
     rho = (options%eta_very_successful + options%eta_too_successful) / 2
     call update_trust_region_radius(rho,options,status,work)
     if ( work%Delta <= 100.0_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp
     
     ! check if rho stays the same because too successful...
     rho = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( abs(work%Delta - 100.0_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp

     rho = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( work%Delta >= 100.0_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp

     rho = options%eta_successful - 10.0_wp
     call update_trust_region_radius(rho,options,status,work)
     if ( work%Delta >= 100.0_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        fails = fails + 1
     end if
     work%Delta = 100.0_wp
     
     ! again...NaN test should go here!!!

     !Finally, check the error cases...
     
     options%tr_update_strategy = 18
     call update_trust_region_radius(rho,options,status,work)
     if ( status%status .ne. NLLS_ERROR_BAD_TR_STRATEGY ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Error returned is = ', status%status, ', expected ',NLLS_ERROR_BAD_TR_STRATEGY
        fails = fails + 1
     end if
     status%status = 0
     work%Delta = 100.0_wp


     
     call reset_default_options(options)
     
   end subroutine update_trust_region_radius_tests

   subroutine test_convergence_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp) :: normF, normJf, normF0, normJF0, normd
     type( nlls_inform ) :: status

     fails = 0
     
     ! hit the case where f = 0 
     normF = 0.0_wp
     normJf = 0.0_wp
     normF0 = 1.0_wp
     normJF0 = 1.0_wp
     normd = 0.0_wp
     call test_convergence(normF,normJf,normF0,normJF0,normd, options, status)
     if ( status%convergence_normf .ne. 1) then
        write(*,*) 'Error in test_convergence test :: expected status%convergence_normf = 1'
        write(*,*) 'got status%convergence_normf = ', status%convergence_normf
        fails = fails + 1
     end if

     call reset_default_options(options)
     
   end subroutine test_convergence_tests

   subroutine mult_J_tests(options,fails)     

     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: J(:), x(:), Jx(:)
     integer :: n, m 

     fails = 0 
     
     n = 2
     m = 4
     allocate(J(n*m), x(m), Jx(n))
     
     x = 1.0_wp
     J = [ 1.0 , 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 ]
     call mult_J(J,m,n,x,Jx)
     if ( norm2( Jx - [16.0, 20.0 ] ) > 1e-12) then
        write(*,*) 'error :: mult_J test failed'
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine mult_J_tests

   subroutine mult_Jt_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: J(:), x(:), Jtx(:)
     integer :: n, m 
     
     fails = 0

     n = 2
     m = 4
     allocate( J(n*m) ) 
     allocate( x(m) ) 
     allocate( Jtx(n) )
     
     x = 1.0_wp
     J = [ 1.0 , 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 ]
     call mult_Jt(J,n,m,x,Jtx)
     if ( norm2( Jtx - [10.0, 26.0 ] ) > 1e-12) then
        write(*,*) 'error :: mult_Jt test failed'
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine mult_Jt_tests


   

   subroutine switch_to_quasi_newton_tests(options,fails)
     
     type( nlls_options), intent(inout) :: options
     integer, intent(out) :: fails

     integer  :: n, m
     TYPE( nlls_inform ) :: inform
     type( nlls_workspace ) :: w
     type( nlls_workspace ), Target :: iw

     fails = 0
     w%iw_ptr => iw
     iw%iw_ptr => iw

     n = 2
     m = 4
     
     options%type_of_method = 2
     options%reg_order = 0.0_wp
     call setup_workspaces(w,n,m,options,inform)
     ! switch_to_quasi_newton expect for w%hf_temp to be defined.
     ! Fill with zeros...
     w%hf_temp(1:n**2) = 0.0_wp
     call switch_to_quasi_newton(w,n,options)
     call remove_workspaces(w,options)

     call reset_default_options(options)
     
   end subroutine switch_to_quasi_newton_tests

   
   subroutine solve_spd_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), b(:), LtL(:,:), x_calc(:), x_true(:)
     integer :: n
     type( nlls_inform ) :: status
     
     fails = 0
    
     n = 2
     allocate(A(n,n), b(n), LtL(n,n), x_calc(n), x_true(n))
     A = reshape([ 4.0, 1.0, 1.0, 2.0 ], shape(A))
     x_true = [1.0, 1.0]
     b = [5.0, 3.0]
     call solve_spd(A,b,LtL,x_calc,n,status)
     if (status%status .ne. 0) then
        write(*,*) 'Error: info = ', status%status, ' returned from solve_spd'
        fails = fails + 1
     else if (norm2(x_calc-x_true) > 1e-12) then
        write(*,*) 'Error: incorrect value returned from solve_spd'
        fails = fails + 1
     end if
          
     call reset_default_options(options)
     
   end subroutine solve_spd_tests

   subroutine solve_general_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), b(:), x_calc(:), x_true(:)
     integer :: n,m
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     
     fails = 0
     work%iw_ptr => iw
     iw%iw_ptr => iw
     
     n = 2
     m = 3
     
     allocate(A(n,n), b(n), x_calc(n), x_true(n))
     
     options%nlls_method = 2
     options%model = 2
     call setup_workspaces(work,n,m,options,status)

     A = reshape([4.0, 1.0, 2.0, 2.0], shape(A))
     x_true = [1.0, 1.0] 
     b = [6.0, 3.0]

     call solve_general(A,b,x_calc,n,status,& 
          work%calculate_step_ws%AINT_tr_ws%solve_general_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: info = ', status%status, ' returned from solve_general'
        fails = fails + 1
        status%status = 0
     else if (norm2(x_true-x_calc) > 1e-12) then
        write(*,*) 'Error: incorrect value returned from solve_general'
        fails = fails + 1
     end if

     A = reshape( [ 0.0, 0.0, 0.0, 0.0 ],shape(A))
     b = [ 6.0, 3.0 ]

     call solve_general(A,b,x_calc,n,status,& 
          work%calculate_step_ws%AINT_tr_ws%solve_general_ws)
     if (status%status .ne. NLLS_ERROR_FROM_EXTERNAL) then
        write(*,*) 'Error: expected error return from solve_general, got info = ', status%status
        fails = fails + 1
     end if

     call nlls_finalize(work,options)
     call reset_default_options(options)
     
   end subroutine solve_general_tests

   subroutine matmult_inner_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), AtA(:,:), AtA_expected(:,:), diff(:)
     integer :: n, m, i

     fails = 0 
     

     n = 2
     m = 3
     allocate(A(m,n), AtA(n,n), AtA_expected(n,n), diff(n))
     A = reshape( [ 1.0, 2.0, 3.0,  &
                    2.0, 4.0, 6.0 ],&
                    shape(A))
     call matmult_inner(A,n,m,AtA)
     AtA_expected = reshape( [ 14.0, 28.0,  &
                               28.0, 56.0 ] &
                               , shape(AtA_expected))
     do i = 1,n
        diff(i) = norm2(AtA(:,i) - AtA_expected(:,i))
     end do
     if (norm2(diff) > 1e-10) then
        write(*,*) 'error :: matmult_inner test failed'
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine matmult_inner_tests

   
   subroutine matmult_outer_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), AAt(:,:), AAt_expected(:,:), diff(:)
     integer :: n, m, i
     
     fails = 0 

     n = 2
     m = 3
     allocate(A(m,n), AAt(m,m), AAt_expected(m,m), diff(m))
     A = reshape( [1.0, 2.0, 3.0,  &
                   2.0, 4.0, 6.0],&
                   shape(A))
     call matmult_outer(A,n,m,AAt)
     AAt_expected = reshape( [  5.0, 10.0, 15.0,  &
                               10.0, 20.0, 30.0, & 
                               15.0, 30.0, 45.0 ] &
                               , shape(AAt_expected))
     do i = 1,m
        diff(i) = norm2(AAt(:,i) - AAt_expected(:,i))
     end do
     if (norm2(diff) > 1e-10) then
        write(*,*) 'error :: matmult_outer test failed'
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine matmult_outer_tests

   subroutine outer_product_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: x(:), xxt(:,:), xxt_exact(:,:), diff(:)
     integer :: n, i 
     
     fails = 0

     n = 4
     allocate(x(n), xxt(n,n), xxt_exact(n,n), diff(n))
     x = [ 1.0, 2.0, 3.0, 4.0 ]
     xxt_exact = reshape( [1.0, 2.0,  3.0,  4.0, &
                           2.0, 4.0,  6.0,  8.0, &
                           3.0, 6.0,  9.0, 12.0, & 
                           4.0, 8.0, 12.0, 16.0], shape(xxt_exact))
     call outer_product(x,n,xxt)
     do i = 1, n
        diff(i) = norm2(xxt(i,:) - xxt_exact(i,:))
     end do
     if (norm2(diff) > 1e-12) then
        write(*,*) 'error :: outer_product test failed'
        fails = fails +1
     end if

     
     call reset_default_options(options)
     
   end subroutine outer_product_tests

   subroutine min_eig_symm_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), ev(:)
     real(wp) :: ew
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     integer :: n, m, i 
     
     fails = 0 
     work%iw_ptr => iw
     iw%iw_ptr => iw
     
     n = 4
     m = 4

     allocate(ev(n), A(n,n))
     
     ! make sure min_eig_symm gets called
     do i = 1, 2
        ! Setup workspace for n = 4
        ! use this for min_eig_symm

        options%nlls_method = 3
        select case (i)
        case (1)
           options%subproblem_eig_fact = .TRUE.
        case (2)
           options%subproblem_eig_fact = .FALSE.
        end select

        ! first, remove previous workspace
        call remove_workspaces(work,options)
        call setup_workspaces(work,n,m,options,status) 
        if (status%status .ne. 0) then
           write(*,*) 'Error: info = ', status%status, ' when setting up workspace'
           fails = fails +1 
        end if
                
        A = reshape( [-5.0, 1.0, 0.0, 0.0, &
                      1.0, -5.0, 0.0, 0.0, &
                      0.0,  0.0, 4.0, 2.0, & 
                      0.0,  0.0, 2.0, 4.0], shape(A))

        call min_eig_symm(A,n,ew,ev,options,status, & 
             work%calculate_step_ws%more_sorensen_ws%min_eig_symm_ws)

        if ( (abs( ew + 6.0 ) > 1e-12).or.(status%status .ne. 0) ) then
           write(*,*) 'error :: min_eig_symm test failed -- wrong eig found'
           fails = fails +1 
        elseif ( norm2(matmul(A,ev) - ew*ev) > 1e-12 ) then
           write(*,*) 'error :: min_eig_symm test failed -- not an eigenvector'
           fails = fails +1 
        end if

        call nlls_finalize(work,options)

        call min_eig_symm(A,n,ew,ev,options,status, & 
             work%calculate_step_ws%more_sorensen_ws%min_eig_symm_ws)
        if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
           write(*,*) 'Error: workspace error not flagged when workspaces not setup'
           fails = fails + 1
        end if
        
     end do

     call reset_default_options(options)
     
   end subroutine min_eig_symm_tests

   subroutine max_eig_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), B(:,:), ev(:), nullevs(:,:), shapenull(:), diff(:)
     real(wp) :: ew
     type( nlls_inform ) :: status
     type( nlls_workspace ) :: work
     type( nlls_workspace ), Target :: iw
     integer :: n, m, i 
     
     fails = 0 
     work%iw_ptr => iw
     iw%iw_ptr => iw
     
     n = 2
     m = 2

     allocate(ev(2*n), A(2*n,2*n), B(2*n,2*n),nullevs(n,2))
     
     options%nlls_method = 2
     call setup_workspaces(work,n,m,options,status)
     
     A = reshape( [ 1.0, 2.0, 3.0,  4.0, &
                    2.0, 4.0, 6.0,  8.0, &
                    3.0, 6.0, 9.0, 12.0, & 
                    4.0, 8.0, 12.0,16.0 ], shape(A))
     B = 0.0_wp
     do i = 1,2*n
        B(i,i) = real(i,wp)
     end do
     !alpha = 1.0_wp
     !x = 0.0_wp
     call max_eig(A,B,2*n,ew,ev,nullevs,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if ( status%status .ne. 0 ) then
        write(*,*) 'error :: max_eig test failed, status = ', status%status
        fails = fails + 1 
     elseif ( (abs( ew - 10.0_wp) > 1e-12) ) then
        write(*,*) 'error :: max_eig test failed, incorrect answer'
        write(*,*) 'expected 10.0, got ', ew
        fails = fails + 1
     end if

     A = 0.0_wp
     A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
     A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
     B = A
     A(1,1) = 1.0_wp; A(2,2) = 1.0_wp
     
     deallocate(nullevs)
     call max_eig(A,B,2*n,ew,ev,nullevs,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if (.not. allocated(nullevs)) then ! check C returned 
        write(*,*) 'error :: hard case of max_eig test failed - C not returned'
        fails = fails + 1
     else
        allocate(shapenull(2))
        shapenull = shape(nullevs)
        if ((shapenull(1) .ne. 2) .or. (shapenull(2) .ne. 2*n)) then
           write(*,*) 'error :: hard case of max_eig test failed - wrong shape C returned'
           write(*,*) 'shapenull(1) = ', shapenull(1), 'shapenull(2) = ', shapenull(2)
           fails = fails + 1
        else
           allocate(diff(n))
           ! Repopulate A (was overwritten by eig routine)
           A = 0.0_wp  
           A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
           A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
           B = A
           A(1,1) = 1.0_wp; A(2,2) = 1.0_wp
           do i = 1, n
              diff(i) = norm2(                        &
                   matmul( A(3:4,3:4),nullevs(1:2,i) )         &
                   - ew * matmul(B(3:4,3:4),nullevs(1:2,i)) & 
                   )
           end do
           if (norm2(diff) > 1e-10) then
              write(*,*) 'error :: hard case of max_eig test failed - wrong vectors returned'
              write(*,*) 'diff = ', diff
              fails = fails + 1
           end if
        end if
     end if
     
     call remove_workspaces(work, options)

     deallocate(ev,A,B)
     
     n = 1
     m = 1
     call setup_workspaces(work, n, m, options, status)

     ! check the error return
     allocate(ev(2*n), A(2*n,2*n), B(2*n,2*n))
     A = 0.0_wp
     B = 0.0_wp
     A(1,2) = 1.0_wp
     A(2,1) = -1.0_wp
     B(1,1) = 1.0_wp
     B(2,2) = 1.0_wp

     call max_eig(A,B,2*n,ew,ev,nullevs,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if (status%status .ne. NLLS_ERROR_AINT_EIG_IMAG) then
        write(*,*) 'error :: all complex part of max_eig test failed'
        fails = fails + 1
     end if
     status%status = 0

     call max_eig(A,B,2*n+1,ew,ev,nullevs, options,status, &
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if ( status%status .ne. NLLS_ERROR_AINT_EIG_ODD ) then
        write(*,*) 'error :: even part of max_eig test failed'
        fails = fails + 1
     end if
     status%status = 0

     call nlls_finalize(work,options)

     ! now let's check for workspace error
     call max_eig(A,B,2*n+1,ew,ev,nullevs, options,status, &
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        fails = fails + 1 
     end if

     
     call reset_default_options(options)
     
   end subroutine max_eig_tests

   subroutine shift_matrix_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: A(:,:), AplusSigma(:,:)
     real(wp) :: sigma
     integer :: n
     
     fails = 0 

     n = 2

     allocate(A(n,n),AplusSigma(n,n))
     A = 1.0_wp
     AplusSigma = 0.0_wp
     sigma = 5.0_wp
     call shift_matrix(A,sigma,AplusSigma,n)
     if ( ( (AplusSigma(1,1)-6.0_wp) > 1e-12) .or. &
          ((AplusSigma(2,2) - 6.0_wp) > 1e-12) ) then
        write(*,*) 'Error: incorrect return from shift_matrix'
        fails = fails + 1
     elseif ( ( (AplusSigma(1,2)-1.0_wp) > 1e-12) .or. & 
          ((AplusSigma(2,1) - 1.0_wp) > 1e-12) ) then
        write(*,*) 'Error: incorrect return from shift_matrix'
        fails = fails + 1
     end if
     
     call reset_default_options(options)
     
   end subroutine shift_matrix_tests
   
   
   subroutine error_message_tests(options,fails)
     
     type( nlls_options ), intent(inout) :: options
     integer, intent(out) :: fails
     
     real(wp), allocatable :: X(:)
     type( user_type ) :: params
     type( nlls_inform ) :: status
     character (len = 80) :: expected_string
     
     integer :: n, m
     

     fails = 0
     n = 2
     m = 67

     allocate(params%x_values(m))
     allocate(params%y_values(m))
     
     call generate_data_example(params)
     
     ! let's check the workspace errors 
     ! first, let's do the main workspace...
     allocate(X(n))
     options%setup_workspaces = .false.
     call nlls_solve(n, m, X,                   &
                    eval_F, eval_J, eval_H, params, &
                    options, status)
     if (status%status .ne. NLLS_ERROR_WORKSPACE_ERROR) then 
        write(*,*) 'Error: workspace error not flagged when workspaces not setup'
        write(*,*) 'status = ', status%status
        fails = fails + 1
     end if

     
     ! nlls_strerror
     status%status = NLLS_ERROR_MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_EVALUATION
     status%external_name = 'nlls_test'
     status%external_return = -1
     write(expected_string,'(a,a,a,i0)') & 
                'Error code from user-supplied subroutine ',trim(status%external_name), & 
                ' passed error = ', status%external_return
     call nlls_strerror(status)
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_UNSUPPORTED_MODEL
     call nlls_strerror(status)
     expected_string = 'Unsupported model passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_FROM_EXTERNAL
     call nlls_strerror(status)
     write(expected_string,'(a,a,a,i0)') & 
                'The external subroutine ',trim(status%external_name), & 
                ' passed error = ', status%external_return
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if
     
     status%status = NLLS_ERROR_UNSUPPORTED_METHOD
     call nlls_strerror(status)
     expected_string = 'Unsupported nlls_method passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_ALLOCATION
     status%bad_alloc = "nlls_test"
     call nlls_strerror(status)
     write(expected_string,'(a,a)') &
                'Bad allocation of memory in ', trim(status%bad_alloc)
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_MAX_TR_REDUCTIONS
     call nlls_strerror(status)
     expected_string = 'The trust region was reduced the maximum number of times'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if
     
     status%status = NLLS_ERROR_MAX_TR_REDUCTIONS
     call nlls_strerror(status)
     expected_string = 'The trust region was reduced the maximum number of times'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_X_NO_PROGRESS
     call nlls_strerror(status)
     expected_string = 'No progress made in X'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_N_GT_M
     call nlls_strerror(status)
     expected_string = 'The problem is overdetermined'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_BAD_TR_STRATEGY
     call nlls_strerror(status)
     expected_string = 'Unsupported tr_update_stategy passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_FIND_BETA
     call nlls_strerror(status)
     expected_string = 'Unable to find suitable scalar in findbeta subroutine'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_BAD_SCALING
     call nlls_strerror(status)
     expected_string = 'Unsupported value of scale passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_WORKSPACE_ERROR
     call nlls_strerror(status)
     expected_string =  'Error accessing pre-allocated workspace'
     if (status%error_message .ne. expected_string) then 
        write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
             status%status
        fails = fails + 1
     end if

     status%status = NLLS_ERROR_UNSUPPORTED_TYPE_METHOD
     call nlls_strerror(status)
     expected_string =  'Unsupported value of type_of_method passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_DOGLEG_MODEL
     call nlls_strerror(status)
     expected_string = 'Model not supported in dogleg (nlls_method=1)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_AINT_EIG_IMAG
     call nlls_strerror(status)
     expected_string = 'All eigenvalues are imaginary (nlls_method=2)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_AINT_EIG_ODD
     call nlls_strerror(status)
     expected_string = 'Odd matrix sent to max_eig subroutine (nlls_method=2)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_MS_MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum iterations reached in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_MS_TOO_MANY_SHIFTS
     call nlls_strerror(status)
     expected_string = 'Too many shifts taken in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if


     status%status = NLLS_ERROR_MS_NO_PROGRESS
     call nlls_strerror(status)
     expected_string = 'No progress being made in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     status%status = NLLS_ERROR_NO_SECOND_DERIVATIVES
     call nlls_strerror(status)
     expected_string = 'Exact second derivatives needed for tensor model'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if
     
     


     status%status = -2355
     call nlls_strerror(status)
     expected_string = 'Unknown error number'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         fails = fails + 1
     end if

     call reset_default_options(options)
     
   end subroutine error_message_tests


   
 end module example_module
