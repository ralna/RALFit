program nlls_test
  
! Test deck for nlls_module

  use ral_nlls_double
  use ral_nlls_internal
  use example_module

  implicit none


  integer, parameter :: wp = kind(1.0d0)
  type( NLLS_inform )  :: status
  type( NLLS_options ) :: options
  type( user_type ), target :: params
  real(wp), allocatable :: v(:),w(:),x(:),y(:),z(:)
  real(wp), allocatable :: A(:,:), B(:,:), C(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha, beta, gamma, delta
  integer :: m, n, i, no_errors_helpers, no_errors_main, info
  integer :: nlls_method, model, tr_update
  logical :: test_all, test_subs
  character (len = 80) :: expected_string

  integer :: number_of_models
  integer, allocatable :: model_to_test(:)

  type( NLLS_workspace ) :: work

  options%error = 18
  options%out   = 17 
  open(unit = options%out, file="nlls_test.out")
  open(unit = options%error, file="nlls_test_error.out")
  
  test_all = .true.
  test_subs = .true.

  no_errors_main = 0

  if (test_all) then
  !!!!!!!!!!!!!!!!!!!!!!!!
  !! Test the main file !!
  !!!!!!!!!!!!!!!!!!!!!!!!
     write(*,*) '==========================='
     write(*,*) '=--Testing the main file--='
     write(*,*) '==========================='

     n = 2

     m = 67
          
     number_of_models = 4
     allocate(model_to_test(number_of_models))
     model_to_test = (/ 0, 1, 2, 3 /)

     allocate( x(n) )

     ! Get params for the function evaluations
     allocate(params%x_values(m))
     allocate(params%y_values(m))
     
     call generate_data_example(params%x_values,params%y_values,m)
     
     
     do tr_update = 1,2
        do nlls_method = 1,4
           do model = 1,number_of_models
     
              X(1) = 1.0 
              X(2) = 2.0

              options%print_level = 3
              options%nlls_method = nlls_method
              options%tr_update_strategy = tr_update
              options%model = model_to_test(model)
              
              call nlls_solve(n, m, X,                         &
                   eval_F, eval_J, eval_H, params,  &
                   options, status )
              if (( nlls_method == 1).and.( options%model > 1)) then
                 if ( status%status .ne. ERROR%DOGLEG_MODEL ) then
                    write(*,*) 'incorrect error return from nlls_solve:'
                    write(*,*) 'NLLS_METHOD = ', nlls_method
                    write(*,*) 'MODEL = ', options%model
                    no_errors_main = no_errors_main + 1
                 end if
              else if ( options%model == 0 ) then
                 if ( status%status .ne. -3 ) then
                    write(*,*) 'incorrect error return from nlls_solve:'
                    write(*,*) 'NLLS_METHOD = ', nlls_method
                    write(*,*) 'MODEL = ', options%model
                    no_errors_main = no_errors_main + 1
                 end if
              else if ( status%status .ne. 0 ) then
                 write(*,*) 'nlls_solve failed to converge:'
                 write(*,*) 'NLLS_METHOD = ', nlls_method
                 write(*,*) 'MODEL = ', options%model
                 write(*,*) 'TR_UPDATE = ', tr_update
                 write(*,*) 'info%status = ', status%status
                 no_errors_main = no_errors_main + 1
              end if

           end do
        end do
     end do
     
     ! Let's get to maxits
     options%maxit = 5
     options%model = 1
     options%nlls_method = 1
      X(1) = 1.0 
      X(2) = 2.0
     call nlls_solve(n, m, X,                         &
                   eval_F, eval_J, eval_H, params,  &
                   options, status )
     if ( status%status .ne. ERROR%MAXITS) then
        write(*,*) 'Error: incorrect error return when maxits expected to be reached'
        
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     options%maxit = 100

     ! Let's get save the arrays
     options%model = 1
     options%nlls_method = 1
     options%output_progress_vectors = .true.
      X(1) = 1.0 
      X(2) = 2.0
     call nlls_solve(n, m, X,                         &
                   eval_F, eval_J, eval_H, params,  &
                   options, status )
     if ( status%status .ne. 0) then
        write(*,*) 'Error: did not converge when output_progress_vectors = true'
        no_errors_main = no_errors_main + 1
        status%status = 0
     end if
     options%output_progress_vectors = .false.

     ! Let's use a relative tr radius
     options%nlls_method = 4
     options%model = 3
     options%exact_second_derivatives = .true.
     options%relative_tr_radius = 1
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     options%relative_tr_radius = 0

     ! Let's use a weighted least squares method
     options%nlls_method = 4
     options%model = 3
     allocate(w(m))
     do i = 1, m
        w(i) = 2.0
     end do
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status, w )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge (weighted):'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     deallocate(w)

     ! Let's do one run with non-exact second derivatives 
     options%nlls_method = 4
     options%model = 3
     options%tr_update_strategy = 1
     options%exact_second_derivatives = .false.
     X = 0.0_wp
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if

     ! And the same with model = 2
     options%nlls_method = 4
     options%model = 2
     options%tr_update_strategy = 1
     options%exact_second_derivatives = .false.
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     
     ! Let's get the method to switch to second-order
     options%nlls_method = 4
     options%model = 3
     options%exact_second_derivatives = .true.
     options%relative_tr_radius = 1
     options%stop_g_absolute = 1e-14
     options%stop_g_relative = 1e-14
     options%print_level = 3
     options%hybrid_tol = 1.0_wp
     params%y_values = params%y_values + 15.0_wp
     X(1) = 7.0; X(2) = -5.0
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     options%relative_tr_radius = 0

     ! Let's get one where ||f|| = 0 
     options%nlls_method = 4
     options%model = 2
     options%relative_tr_radius = 1
     options%stop_g_absolute = 1e-4
     options%stop_g_relative = 0.0_wp!1e-4
     options%print_level = 3
     params%y_values = exp( 0.3_wp * params%x_values + 0.2_wp)
     X(1) = 0.3; X(2) = 0.1
     call nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     options%relative_tr_radius = 0
     
     ! now let's check errors on the parameters passed to the routine...
     
     options%print_level = 3

     ! test n > m
     n = 100
     m = 3
     call  nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if (status%status .ne. ERROR%N_GT_M) then
        write(*,*) 'Error: wrong error return, n > m'
        no_errors_main = no_errors_main + 1
     end if
     n = 2
     m = 67
     
    ! test for unsupported method
     options%nlls_method = 3125
     call nlls_solve(n, m, X,                   &
                    eval_F, eval_J, eval_H, params, &
                    options, status)
     if ( status%status .ne. ERROR%UNSUPPORTED_METHOD ) then 
        write(*,*) 'Error: unsupported method passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     options%nlls_method = 4

     ! test for unsupported tr strategy
     options%tr_update_strategy = 323
     call nlls_solve(n, m, X,                   &
                    eval_F, eval_J, eval_H, params, &
                    options, status)
     if ( status%status .ne. ERROR%BAD_TR_STRATEGY ) then 
        write(*,*) 'Error: unsupported TR strategy passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     options%tr_update_strategy = 1
     
     if (no_errors_main == 0) then
        write(*,*) '*** All (main) tests passed successfully! ***'
     else
        write(*,*) 'There were ', no_errors_main,' errors'
     end if

     deallocate(x,params%x_values,params%y_values)
     
     
  end if

  deallocate(model_to_test)

  no_errors_helpers = 0
  
  if ( test_subs ) then 

     !###############################!
     !###############################!
     !! Test the helper subroutines !!
     !###############################!
     !###############################!

     write(*,*) '============================='
     write(*,*) '=--Testing the subroutines--='
     write(*,*) '============================='
     

     options%print_level = 3
     options%scale = 0

     !! calculate step...
     ! not needed -- fully tested elsewhere....

     !! dogleg 
     options%nlls_method = 1
     options%model = 5
     n = 2
     m = 3
     allocate(w(m*n), x(n*n), y(m), z(n), v(n))
     call setup_workspaces(work,n,m,options,status) 
     ! w <-- J
     ! x <-- hf
     ! y <-- f
     ! z <-- g 
     ! v <-- d
     alpha = 10.0_wp
     
     ! first, hit the 'method not supported' error
     options%model = 27
     w = 1.0_wp
     x = 0.0_wp
     y = 1.0_wp
     z = 1.0_wp
     call dogleg(w,y,x,z,n,m,alpha,v,beta,options,status,work%calculate_step_ws%dogleg_ws)
     if (status%status .ne. ERROR%DOGLEG_MODEL) then
        write(*,*) 'Error: unsupported model allowed in dogleg'
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     options%model = 1

     w = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
     x = 0.0_wp
     y = 1.0_wp
     z = 1.0_wp
     ! now, get ||d_gn|| <= Delta
     alpha = 6.0_wp
     call dogleg(w,y,x,z,n,m,alpha,v,beta,options,status,work%calculate_step_ws%dogleg_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in dogleg'
        no_errors_helpers = no_errors_helpers + 1
     end if
     ! now set delta so that || alpha * d_sd || >= Delta
     alpha = 0.5_wp
     call dogleg(w,y,x,z,n,m,alpha,v,beta,options,status,work%calculate_step_ws%dogleg_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in dogleg'
        no_errors_helpers = no_errors_helpers + 1
     end if
     ! now get the guys in the middle...
     alpha = 2.5_wp
     call dogleg(w,y,x,z,n,m,alpha,v,beta,options,status,work%calculate_step_ws%dogleg_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in dogleg'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     end if
     
     deallocate(x,y,z,v,w)
     call nlls_finalize(work,options)
     
     options%scale = 4
     options%nlls_method = 3
     n = 2
     m = 3 
     allocate(w(m*n),A(n,n),y(n))
     call setup_workspaces(work,n,m,options,status)
     
     w = 0.0_wp
     w(1) = 1e15
     w(4) = 1e-15
     A(1,1) = 1.0_wp
     A(2,1) = 0.0_wp
     A(1,2) = 0.0_wp
     A(2,2) = 1.0_wp

     !** scale = 1 **
     options%scale= 1     
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error in apply_scaling when scale = 1'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if

!!$     !** scale = 2 **
!!$     options%scale = 2
!!$     call apply_scaling(w,n,m,A,y,& 
!!$          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
!!$          options,status)
!!$     if (status%status .ne. 0 ) then
!!$        write(*,*) 'Error: unexpected error in apply_scaling when scale = 2'
!!$        write(*,*) 'status = ', status%status,' returned.'
!!$        no_errors_helpers = no_errors_helpers + 1
!!$        status%status = 0 
!!$     end if

     !** scale = 2 **
     options%scale = 2
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error in apply_scaling when scale = 2'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if

!!$     !** scale = 3 **
!!$     options%scale = 3
!!$     call apply_scaling(w,n,m,A,y,& 
!!$          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
!!$          options,status)
!!$     if (status%status .ne. 0 ) then
!!$        write(*,*) 'Error: unexpected error in apply_scaling when scale = 3'
!!$        write(*,*) 'status = ', status%status,' returned.'
!!$        no_errors_helpers = no_errors_helpers + 1
!!$        status%status = 0 
!!$     end if

     !** scale undefined
     options%scale = 786
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. ERROR%BAD_SCALING ) then
        write(*,*) 'Error: expected error in apply_scaling when passing undefined scaling'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if
     status%status = 0
     
     ! now, let's test the non-default modes
     ! first, set scale_require_increase to T
     options%scale = 1
     options%scale_require_increase = .true.
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if
     options%scale_require_increase = .false.

     ! first, set scale_trim_min to T
     options%scale_trim_min = .true.
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if
     options%scale_trim_min = .false.

     ! first, set scale_trim_max to T
     options%scale_trim_max = .false.
     call apply_scaling(w,n,m,A,y,& 
          work%calculate_step_ws%more_sorensen_ws%apply_scaling_ws, &
          options,status)
     if (status%status .ne. 0 ) then
        write(*,*) 'Error: unexpected error when scale_require_increase = T'
        write(*,*) 'status = ', status%status,' returned.'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0 
     end if
     options%scale_trim_max = .true.

     
     deallocate(w,A,y)
     call nlls_finalize(work,options)
     options%scale = 0 

     !! aint_tr
     ! ** TODO ** 

     !! more_sorensen
     options%nlls_method = 3
     n = 2
     m = 3
     allocate(w(m*n), x(n*n), y(m), z(n))
     call setup_workspaces(work,n,m,options,status) 
     ! w <-- J
     ! x <-- hf
     ! y <-- f
     ! z <-- d 
     alpha = 10.0_wp
     
     ! regular case...
     w = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
     x = 0.0_wp
     y = 1.0_wp
     z = 1.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     end if

     ! non spd matrix, with failure
     options%more_sorensen_shift = -1000.0_wp
     w = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
     x = -10 * (/ 2.0_wp, 1.0_wp, 5.0_wp, 7.0_wp /)
     y = 1.0_wp
     z = 1.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. ERROR%MS_TOO_MANY_SHIFTS) then
        write(*,*) 'Error: MS too many shifts test passed, when fail expected'
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     options%more_sorensen_shift = 1e-13

     ! look for nd /=  Delta with a non-zero shift?
     w = (/ 2.0_wp, 3.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 4.0_wp /)
     x = (/ 1.0_wp, 2.0_wp, 2.0_wp, 1.0_wp /)
     y = 1.0_wp
     z = 1.0_wp
     alpha =  10.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with non-zero shift'
        write(*,*) 'status = ', status%status, ' returned'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     end if

     ! now look for nd =  Delta with a non-zero shift?
     w = (/ 2.0_wp, 3.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 4.0_wp /)
     x = (/ 1.0_wp, 2.0_wp, 2.0_wp, 1.0_wp /)
     y = 1.0_wp
     z = 1.0_wp
     beta = options%more_sorensen_tiny
     options%more_sorensen_tiny = 0.01_wp
     alpha =  0.2055_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with non-zero shift'
        write(*,*) 'status = ', status%status, ' returned'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     end if
     options%more_sorensen_tiny = beta
     ! *todo*

     ! now take nd > Delta
     w = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
     x = 0.0_wp
     y = 1.0_wp
     z = 1.0_wp
     alpha = 3.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: unexpected error in more-sorensen test with nd > Delta'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     end if


     ! get to max_its...
     options%more_sorensen_maxits = 1     
     w = 0.1_wp * (/ 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
     x = 0.0_wp
     y = 1.0_wp
     z = 1.0_wp
     alpha = 3.0_wp
     ! now, get ||d_gn|| <= Delta
     call more_sorensen(w,y,x,n,m,alpha,z,beta,options,status,& 
          work%calculate_step_ws%more_sorensen_ws)
     if (status%status .ne. ERROR%MS_MAXITS) then
        write(*,*) 'Error: Expected maximum iterations error in more_sorensen'
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     options%more_sorensen_maxits = 10
     
     deallocate(x,y,z,w)
     call nlls_finalize(work,options)

     !! solve_galahad
     do i = 1,2
        options%type_of_method = i
        options%nlls_method = 4
        n = 2
        m = 5
        call setup_workspaces(work,n,m,options,status) 

        allocate(w(n))
        allocate(x(m*n))
        allocate(y(m))
        allocate(z(n*n))
        ! x -> J, y-> f, x -> hf, w-> d
        x = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp /)
        y = (/ 1.2_wp, 3.1_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
        z = 1.0_wp

        alpha = 0.02_wp
        
        work%calculate_step_ws%solve_galahad_ws%reg_order = 2.0_wp
        call solve_galahad(x,y,z,n,m,alpha,w,beta,& 
             options,status, &
             work%calculate_step_ws%solve_galahad_ws )
        
        if ( status%status .ne. 0 ) then
           select case (i)
              case(1)
                 write(*,*) 'DTRS test failed, status = ', status%status
              case(2)
                 write(*,*) 'DRQS test failed, status = ', status%status
              end select
              no_errors_helpers = no_errors_helpers + 1
        end if
     
        if (i == 1) then
           ! check result lies within the trust region
           if ( abs(dot_product(w,w) - alpha**2) > 1e-12 ) then
              select case (i)
              case (1)
                 write(*,*) 'dtrs failed'
              case (2)
                 write(*,*) 'drqs failed'
              end select
              write(*,*) 'Delta = ', alpha, '||d|| = ', dot_product(w,w)
              no_errors_helpers = no_errors_helpers + 1
           end if
        end if

        ! Flag an error from dtrs...
        x = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp /)
        y = (/ 1.2_wp, 3.1_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
        z = 1.0_wp

        alpha = -100.0_wp
        
        call solve_galahad(x,y,z,n,m,alpha,w,beta,& 
             options,status,& 
             work%calculate_step_ws%solve_galahad_ws)

        if ( status%status .ne. ERROR%FROM_EXTERNAL ) then
           select case (i)
           case (1)
              write(*,*) 'DTRS test failed, expected status = ', ERROR%FROM_EXTERNAL
           case (2)
              write(*,*) 'DRQS test failed, expected status = ', ERROR%FROM_EXTERNAL
           end select
           write(*,*) ' but got status = ', status%status
           no_errors_helpers = no_errors_helpers + 1
        end if
        status%status = 0

     
        deallocate(x,y,z,w)
        call nlls_finalize(work,options)
     end do
     !! solve_LLS 
     options%nlls_method = 1 ! dogleg
     call setup_workspaces(work,n,m,options,status) 
     

     n = 2 
     m = 5
     allocate(x(n*m), y(n), w(m), z(m))
     ! x<--J
     ! z<--f
     ! y<--sol
     ! w<--J*sol
     x = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, & 
            6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp /)
     z = (/ 7.0_wp, 9.0_wp, 11.0_wp, 13.0_wp, 15.0_wp /)
     
     call solve_LLS(x,z,n,m,y,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if ( status%status .ne. 0 ) then 
        write(*,*) 'solve_LLS test failed: wrong error message returned'
        write(*,*) 'status = ', status%status
        no_errors_helpers = no_errors_helpers+1
     end if
     call mult_J(x,n,m,y,w)
     alpha = norm2(w + z)
     if ( alpha > 1e-12 ) then
        ! wrong answer, as data chosen to fit
        write(*,*) 'solve_LLS test failed: wrong solution returned'
        write(*,*) '||Jx - f|| = ', alpha
        no_errors_helpers = no_errors_helpers+1
     end if
     
     ! finally, let's flag an error....
     deallocate(w,x,y,z)
     call nlls_finalize(work, options)
     
     n = 100 
     m = 20
     allocate(x(n*m), y(m), z(n))     
     call setup_workspaces(work,n,m,options,status) 

     x = 1.0_wp
     z = 1.0_wp
     call solve_LLS(x,z,n,m,y,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if ( status%status .ne. ERROR%FROM_EXTERNAL ) then 
        write(*,*) 'solve_LLS test failed: wrong error message returned'
        write(*,*) 'status = ', status%status
        no_errors_helpers = no_errors_helpers+1
     end if
     status%status = 0
     
     deallocate(x,y,z)
     call nlls_finalize(work, options)
     options%nlls_method = 9 ! back to hybrid
     
     !------------!
     !! findbeta !!
     !------------!
     n = 3
     allocate(x(n),y(n),z(n))

     x = (/ 1.0, 2.0, 3.0 /) 
     y = (/ 2.0, 1.0, 1.0 /)

     call findbeta(x,y,10.0_wp,alpha,status)

     if (status%status .ne. 0) then
        write(*,*) 'error -- findbeta did not work: info /= 0'
        no_errors_helpers = no_errors_helpers + 1
     else if ( ( norm2( x + alpha * y ) - 10.0_wp ) > 1e-12 ) then
        write(*,*) 'error -- findbeta did not work'
        write(*,*) '|| x + beta y|| = ', norm2( (x + alpha * y)-10.0_wp)
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     deallocate(x,y,z)
     
     n = 2
     allocate(x(n),y(n),z(n))
     
     x = 1e8_wp
     y = 1.0_wp
     alpha = 1e6
     beta = 0.0_wp

     call findbeta(x,y,beta,gamma,status)

     if (status%status .ne. ERROR%FIND_BETA) then
        write(*,*) 'Expected an error from findbeta: info =', status%status
        write(*,*) 'beta returned = ', gamma
        write(*,*) '|| x + beta y|| = ', norm2( (alpha * x + gamma * y)-beta)
        no_errors_helpers = no_errors_helpers + 1
     end if

     deallocate(x,y,z)

     !------------------!
     !! evaluate_model !!
     !------------------!

     !! todo
     
     !-----------------!
     !! calculate_rho !!
     !-----------------!

     alpha = 2.0_wp ! normf
     beta =  1.0_wp ! normfnew
     gamma = 1.5_wp ! md
     call calculate_rho(alpha, beta, gamma, delta,options)
     if ( abs(delta - 3.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 3.0, got ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     ! now, let's check one is returned if alpha = beta
     beta = 2.0_wp
     call calculate_rho(alpha,beta,gamma, delta, options)
     if (abs(delta - 1.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 1.0, got ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     beta = 1.0_wp

     ! finally, check that 1 is returned if denominator = 0
     gamma = 2.0_wp
     call calculate_rho(alpha,beta,gamma, delta, options)
     if (abs(delta - 1.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 1.0, got ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     !! Apply second order info
     ! todo

     !------------------------------!
     !! update_trust_region_radius !!
     !------------------------------!
     call setup_workspaces(work,2,2,options,status) 
     work%Delta = 100.0_wp ! Delta
     work%tr_nu = 2.0_wp ! nu
     work%tr_p = 3 ! p
     ! alpha = rho

     options%tr_update_strategy = 1
     ! let's go through the options
     
     options%eta_success_but_reduce = 0.25_wp
     options%eta_very_successful = 0.75_wp
     options%eta_too_successful = 2.0_wp

     ! check if rho reduced...
     alpha = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( work%Delta >= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp
     
     ! check if rho stays the same...
     alpha = (options%eta_success_but_reduce + options%eta_very_successful) / 2
     call update_trust_region_radius(alpha,options,status,work)
     if ( abs(work%Delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: Delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp

     ! check if rho increases...
     alpha = (options%eta_very_successful + options%eta_too_successful) / 2
     work%normd = 100_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( work%Delta <= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp

     
     ! check if rho stays the same because too successful...
     alpha = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( abs(work%Delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp

     ! now check for NaNs...HOW to do this in a non-compiler dependent way!?!?

     !! now, let's check the other option....
     options%tr_update_strategy = 2
     
     ! check if rho increases...
     alpha = (options%eta_very_successful + options%eta_too_successful) / 2
     call update_trust_region_radius(alpha,options,status,work)
     if ( work%Delta <= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp
     
     ! check if rho stays the same because too successful...
     alpha = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( abs(work%Delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp

     alpha = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( work%Delta >= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp

     alpha = options%eta_successful - 10.0_wp
     call update_trust_region_radius(alpha,options,status,work)
     if ( work%Delta >= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', work%Delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     work%Delta = 100_wp
     
     ! again...NaN test should go here!!!

     !Finally, check the error cases...
     
     options%tr_update_strategy = 18
     call update_trust_region_radius(alpha,options,status,work)
     if ( status%status .ne. ERROR%BAD_TR_STRATEGY ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Error returned is = ', status%status, ', expected ',ERROR%BAD_TR_STRATEGY
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     work%Delta = 100_wp

     !! test_convergence
     
     ! hit the case where f = 0 
     call test_convergence(0.0_wp,1.0_wp, 1.0_wp, 1.0_wp, options, status)
     if ( status%convergence_normf .ne. 1) then
        write(*,*) 'Error in test_convergence test :: expected status%convergence_normf = 1'
        write(*,*) 'got status%convergence_normf = ', status%convergence_normf
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     
     !----------!
     !! mult_J !!
     !----------!

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

     !-----------!
     !! mult_Jt !!
     !-----------!

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

!!!!!!
     ! Setup workspace for n = 2
     ! use this for max_eig, solve_spd
     options%nlls_method = 2
     call setup_workspaces(work,2,2,options,status) 
!!!!!!

     !-------------!
     !! solve_spd !!
     !-------------!

     n = 2
     allocate(A(n,n),x(2),y(2),z(2),B(n,n))
     A = reshape((/ 4.0, 1.0, 1.0, 2.0 /),shape(A))
     z = (/ 1.0, 1.0 /)
     y = (/ 5.0, 3.0 /)

     call solve_spd(A,y,B,x,n,status)
     if (status%status .ne. 0) then
        write(*,*) 'Error: info = ', status%status, ' returned from solve_spd'
        no_errors_helpers = no_errors_helpers + 1
     else if (norm2(x-z) > 1e-12) then
        write(*,*) 'Error: incorrect value returned from solve_spd'
        no_errors_helpers = no_errors_helpers + 1
     end if

     deallocate(A,B,x,y,z)

     !-----------------!
     !! solve_general !!
     !-----------------!
     
     n = 2
     m =2
     options%nlls_method = 2
     options%model = 2

     call setup_workspaces(work,2,2,options,status) 
     allocate(A(n,n),x(2),y(2),z(2))

     A = reshape((/ 4.0, 1.0, 2.0, 2.0 /),shape(A))
     z = (/ 1.0, 1.0 /)
     y = (/ 6.0, 3.0 /)

     call solve_general(A,y,x,n,status,& 
          work%calculate_step_ws%AINT_tr_ws%solve_general_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: info = ', info, ' returned from solve_general'
        no_errors_helpers = no_errors_helpers + 1
        status%status = 0
     else if (norm2(x-z) > 1e-12) then
        write(*,*) 'Error: incorrect value returned from solve_general'
        no_errors_helpers = no_errors_helpers + 1
     end if

     A = reshape((/ 0.0, 0.0, 0.0, 0.0 /),shape(A))
     z = (/ 1.0, 1.0 /)
     y = (/ 6.0, 3.0 /)

     call solve_general(A,y,x,n,status,& 
          work%calculate_step_ws%AINT_tr_ws%solve_general_ws)
     if (status%status .ne. ERROR%FROM_EXTERNAL) then
        write(*,*) 'Error: expected error return from solve_general, got info = ', info
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     
     deallocate(A,x,y,z)
     call nlls_finalize(work,options)
     
     !---------------!
     !! matrix_norm !!
     !---------------!

     ! todo

     !-----------------!
     !! matmult_inner !!
     !-----------------!

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


     !-----------------!
     !! matmult_outer !!
     !-----------------!

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

     !-----------------!
     !! outer_product !!
     !-----------------!

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

     !! All_eig_symm
     ! todo
     
     !----------------!
     !! min_eig_symm !!
     !----------------!

     n = 4
     m = 4
     allocate(x(n),A(n,n))

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
        call setup_workspaces(work,n,m,options,status) 
                
        A = reshape( (/-5.0,  1.0, 0.0, 0.0, &
          1.0, -5.0, 0.0, 0.0, &
          0.0,  0.0, 4.0, 2.0, & 
          0.0,  0.0, 2.0, 4.0/), shape(A))

        call min_eig_symm(A,n,alpha,x,options,status, & 
             work%calculate_step_ws%more_sorensen_ws%min_eig_symm_ws)

        if ( (abs( alpha + 6.0 ) > 1e-12).or.(status%status .ne. 0) ) then
           write(*,*) 'error :: min_eig_symm test failed -- wrong eig found'
           no_errors_helpers = no_errors_helpers + 1 
        elseif ( norm2(matmul(A,x) - alpha*x) > 1e-12 ) then
           write(*,*) 'error :: min_eig_symm test failed -- not an eigenvector'
           no_errors_helpers = no_errors_helpers + 1
        end if

        call nlls_finalize(work,options)
        options%nlls_method = 2 ! revert...

     end do

     deallocate(A,x)
     
!!$     n = 3
!!$     m = 3
!!$     allocate(x(n),A(n,n))
!!$     call setup_workspaces(work,n,m,options,info) 
!!$     options%subproblem_eig_fact = .TRUE.
!!$     
!!$     A = reshape( (/ 1674.456299, -874.579834,  -799.876465,
!!$                     -874.579834,  1799.875854, -925.296021,
!!$                     -799.876465,  -925.296021, 1725.172485/), 
!!$                     shape(A))
!!$
!!$     options%nlls_method = 3
!!$
!!$     call min_eig_symm(A,n,alpha,x,options,status, & 
!!$             work%calculate_step_ws%more_sorensen_ws%min_eig_symm_ws)
!!$     
!!$     call nlls_finalize(work,options)
!!$     deallocate(A,x)    


     !-----------!
     !! max_eig !!
     !-----------!
     n = 2
     m = 2
     ! make sure max_eig gets called

     allocate(x(2*n),A(2*n,2*n), B(2*n,2*n))
     call setup_workspaces(work,n,m,options,status) 
     
     A = reshape( (/1.0, 2.0, 3.0, 4.0, &
          2.0, 4.0, 6.0, 8.0, &
          3.0, 6.0, 9.0, 12.0, & 
          4.0, 8.0, 12.0, 16.0/), shape(A))
     B = 0.0_wp
     do i = 1,2*n
        B(i,i) = real(i,wp)
     end do
     alpha = 1.0_wp
     x = 0.0_wp
     call max_eig(A,B,2*n,alpha,x,C,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if ( status%status .ne. 0 ) then
        write(*,*) 'error :: max_eig test failed, status = ', status%status
        no_errors_helpers = no_errors_helpers + 1 
     elseif ( (abs( alpha - 10.0_wp) > 1e-12) ) then
        write(*,*) 'error :: max_eig test failed, incorrect answer'
        write(*,*) 'expected 10.0, got ', alpha
        no_errors_helpers = no_errors_helpers + 1 
     end if

     deallocate(A,B,x)
     if(allocated(C)) deallocate(C)
     call nlls_finalize(work,options)

     ! check the 'hard' case...
     n = 2
     m = 2
     allocate(x(2*n),A(2*n,2*n), B(2*n,2*n))
     call setup_workspaces(work,n,m,options,status) 

     A = 0.0_wp  
     A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
     A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
     B = A
     A(1,1) = 1.0_wp; A(2,2) = 1.0_wp

     call max_eig(A,B,2*n,alpha,x,C,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if (.not. allocated(C)) then ! check C returned 
        write(*,*) 'error :: hard case of max_eig test failed - C not returned'
        no_errors_helpers = no_errors_helpers + 1 
     else
        allocate(y(2))
        y = shape(C)
        if ((y(1) .ne. 2) .or. (y(2) .ne. 2*n)) then
           write(*,*) 'error :: hard case of max_eig test failed - wrong shape C returned'
           write(*,*) 'y(1) = ', y(1), 'y(2) = ', y(2)
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

     if(allocated(A)) deallocate(A)
     if(allocated(B)) deallocate(B)
     if(allocated(C)) deallocate(C)
     if(allocated(x)) deallocate(x)
     if(allocated(y)) deallocate(y)
     if (allocated(results)) deallocate(results)
     call nlls_finalize(work,options)

     
     
     call setup_workspaces(work,1,1,options,status)  !todo: deallocation routine
     ! check the error return
     n = 2
     allocate(x(n), A(n,n), B(n,n))
     A = 0.0_wp
     B = 0.0_wp
     A(1,2) = 1.0_wp
     A(2,1) = -1.0_wp
     B(1,1) = 1.0_wp
     B(2,2) = 1.0_wp

     call max_eig(A,B,n,alpha,x,C,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if (status%status .ne. ERROR%AINT_EIG_IMAG) then
        write(*,*) 'error :: all complex part of max_eig test failed'
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0

     call max_eig(A,B,n+1,alpha,x,C, options,status, &
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)
     if ( status%status .ne. ERROR%AINT_EIG_ODD ) then
        write(*,*) 'error :: even part of max_eig test failed'
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0

     deallocate(A,B,x)
     call nlls_finalize(work,options)

     !! shift_matrix 
     n = 2
     allocate(A(2,2),B(2,2))
     A = 1.0_wp
     B = 0.0_wp
     alpha = 5.0_wp
     call shift_matrix(A,alpha,B,n)
     if ( ( (B(1,1)-6.0_wp) > 1e-12) .or. ((B(2,2) - 6.0_wp) > 1e-12) ) then
        write(*,*) 'Error: incorrect return from shift_matrix'
        no_errors_helpers = no_errors_helpers + 1
     elseif ( ( (B(1,2)-1.0_wp) > 1e-12) .or. ((B(2,1) - 1.0_wp) > 1e-12) ) then
        write(*,*) 'Error: incorrect return from shift_matrix'
        no_errors_helpers = no_errors_helpers + 1
     end if
     deallocate(A,B)
     
     !! get_svd_J 
     ! Todo

     !! let's make sure output_progress_vectors gets hit
     options%output_progress_vectors = .true.

     n = 2
     m = 3
     call setup_workspaces(work,n,m,options,status)    
     call nlls_finalize(work,options)
     
     ! nlls_strerror
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%EVALUATION
     status%external_name = 'nlls_test'
     status%external_return = -1
     write(expected_string,'(a,a,a,i0)') & 
                'Error code from user-supplied subroutine ',trim(status%external_name), & 
                ' passed error = ', status%external_return
     call nlls_strerror(status)
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%UNSUPPORTED_MODEL
     call nlls_strerror(status)
     expected_string = 'Unsupported model passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%FROM_EXTERNAL
     call nlls_strerror(status)
     write(expected_string,'(a,a,a,i0)') & 
                'The external subroutine ',trim(status%external_name), & 
                ' passed error = ', status%external_return
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = ERROR%MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum number of iterations reached'
     if (status%error_message .ne. 'Maximum number of iterations reached') then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     
     status%status = ERROR%UNSUPPORTED_METHOD
     call nlls_strerror(status)
     expected_string = 'Unsupported nlls_method passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%ALLOCATION
     status%bad_alloc = "nlls_test"
     call nlls_strerror(status)
     write(expected_string,'(a,a)') &
                'Bad allocation of memory in ', trim(status%bad_alloc)
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     status%status = ERROR%MAX_TR_REDUCTIONS
     call nlls_strerror(status)
     expected_string = 'The trust region was reduced the maximum number of times'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if
     
     status%status = ERROR%MAX_TR_REDUCTIONS
     call nlls_strerror(status)
     expected_string = 'The trust region was reduced the maximum number of times'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%X_NO_PROGRESS
     call nlls_strerror(status)
     expected_string = 'No progress made in X'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%N_GT_M
     call nlls_strerror(status)
     expected_string = 'The problem is overdetermined'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%BAD_TR_STRATEGY
     call nlls_strerror(status)
     expected_string = 'Unsupported tr_update_stategy passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%FIND_BETA
     call nlls_strerror(status)
     expected_string = 'Unable to find suitable scalar in findbeta subroutine'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%BAD_SCALING
     call nlls_strerror(status)
     expected_string = 'Unsupported value of scale passed in options'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%DOGLEG_MODEL
     call nlls_strerror(status)
     expected_string = 'Model not supported in dogleg (nlls_method=1)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%AINT_EIG_IMAG
     call nlls_strerror(status)
     expected_string = 'All eigenvalues are imaginary (nlls_method=2)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%AINT_EIG_ODD
     call nlls_strerror(status)
     expected_string = 'Odd matrix sent to max_eig subroutine (nlls_method=2)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%MS_MAXITS
     call nlls_strerror(status)
     expected_string = 'Maximum iterations reached in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%MS_TOO_MANY_SHIFTS
     call nlls_strerror(status)
     expected_string = 'Too many shifts taken in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = ERROR%MS_NO_PROGRESS
     call nlls_strerror(status)
     expected_string = 'No progress being made in more_sorensen (nlls_method=3)'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if


     status%status = -2355
     call nlls_strerror(status)
     expected_string = 'Unknown error number'
     if (status%error_message .ne. expected_string) then 
         write(*,*) 'Error: incorrect string returned from nlls_strerror when status = ', &
              status%status
         no_errors_helpers = no_errors_helpers + 1
     end if

     ! Report back results....

     if (no_errors_helpers == 0) then
        write(*,*) '*** All (helper) tests passed successfully! ***'
     else
        write(*,*) 'There were ', no_errors_helpers,' errors'
     end if

  end if

  
close(unit = 17)
!
!no_errors_helpers = 1
 if (no_errors_helpers + no_errors_main == 0) then
    write(*,*) ' '
    write(*,*) '**************************************'
    write(*,*) '*** All tests passed successfully! ***'
    write(*,*) '**************************************'
    stop 0    ! needed for talking with ctest
 else 
    stop 1    ! needed for talking with ctest
  end if
  


end program nlls_test
