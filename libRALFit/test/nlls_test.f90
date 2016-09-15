program nlls_test
  
! Test deck for nlls_module

! use ral_nlls_double
! use ral_nlls_internal
  use example_module

  implicit none


  type( NLLS_inform )  :: status
  type( NLLS_options ) :: options
  type( user_type ), target :: params
  real(wp), allocatable :: v(:),w(:),x(:),y(:),z(:)
  real(wp), allocatable :: A(:,:), B(:,:), C(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha, beta, gamma, delta
  integer :: m, n, i, no_errors_helpers, no_errors_main, info
  integer :: nlls_method, model, tr_update, inner_method
  logical :: test_all, test_subs
  character (len = 80) :: expected_string
  integer :: fails

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
     
!!$     ! Get params for the function evaluations
!!$     allocate(params%x_values(m))
!!$     allocate(params%y_values(m))
!!$     
!!$     call generate_data_example(params%x_values,params%y_values,m)      
     call generate_data_example(params)
     options%print_level = 3     

     do tr_update = 1,2
        do nlls_method = 1,4
           do model = 1,number_of_models

              call reset_default_options(options)
              options%nlls_method = nlls_method
              options%tr_update_strategy = tr_update
              options%model = model_to_test(model)
              call solve_basic(X,params,options,status)
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
                 write(*,*) status%error_message
                 write(*,*) 'NLLS_METHOD = ', nlls_method
                 write(*,*) 'MODEL = ', options%model
                 write(*,*) 'TR_UPDATE = ', tr_update
                 write(*,*) 'info%status = ', status%status
                 write(*,*) 'scale? ', options%scale
                 no_errors_main = no_errors_main + 1
              end if
           end do
        end do
     end do
     
     ! now, let's test the regularization method
     call reset_default_options(options)
     options%type_of_method = 2 ! regularization
     options%inner_method = 2 
     options%print_level = 1
     options%exact_second_derivatives = .true.
     options%calculate_svd_J = .true.
     
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: regularization'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     
     ! Let's do a test where the regularization weight is non-zero
     call reset_default_options(options)
     options%regularization = 1
     options%regularization_term = 1e-2
     options%regularization_power = 2.0_wp
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: non-zero regularization weight'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! and, the same test with a regularization power of three:
     call reset_default_options(options)
     options%regularization = 1
     options%regularization_term = 1e-2
     options%regularization_power = 3.0_wp
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: regularization_power = 3.0'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     
     ! now let's get regularization with model = 2
     call reset_default_options(options)
     options%regularization = 1
     options%regularization_term = 1e-2
     options%regularization_power = 3.0_wp
     options%model = 2
     options%exact_second_derivatives = .true.
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: regularization_power = 3.0'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     
     ! now test optimal reg_order
     call reset_default_options(options)
     options%type_of_method = 2
     options%model = 1
     options%reg_order = -1.0
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: negative reg_order'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! now, let's do the tensor model...
     call reset_default_options(options)
     options%type_of_method = 2 ! regularization
     options%model = 4 ! hybrid model
     options%exact_second_derivatives = .true.

     do inner_method = 1,2
        options%inner_method = inner_method
        call solve_basic(X,params,options,status)
        if ( status%status .ne. 0 ) then
           write(*,*) 'nlls_solve failed to converge: tensor model'
           write(*,*) 'NLLS_METHOD = ', options%nlls_method
           write(*,*) 'MODEL = ', options%model
           write(*,*) 'info%status = ', status%status
           no_errors_main = no_errors_main + 1
        end if
     end do

     ! now, let's get an error return...
     call reset_default_options(options)
     options%model = 4
     options%exact_second_derivatives = .false.
     call solve_basic(X,params,options,status)
     if ( status%status .ne. ERROR%NO_SECOND_DERIVATIVES ) then
        write(*,*) 'expected error return', ERROR%NO_SECOND_DERIVATIVES,' but'
        write(*,*) 'got ', status%status
        no_errors_main = no_errors_main + 1
     end if


     ! Let's get a subproblem solver error
     call reset_default_options(options)
     options%model = 4
     options%nlls_method = 1
     options%type_of_method = 2
     call solve_basic(X,params,options,status)
     if ( status%status .ne. ERROR%NT_BAD_SUBPROBLEM) then
        write(*,*) 'Error: incorrect error return when wrong subproblem solver selected'
        write(*,*) 'Expected', ERROR%NT_BAD_SUBPROBLEM, ' but got ', status%status
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0 

     ! Let's get to maxits
     call reset_default_options(options)
     options%type_of_method = 1
     options%maxit = 5
     options%model = 1
     options%nlls_method = 1
     call solve_basic(X,params,options,status)
     if ( status%status .ne. ERROR%MAXITS) then
        write(*,*) 'Error: incorrect error return when maxits expected to be reached'
        write(*,*) 'status%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     options%maxit = 100

     ! Let's get save the arrays
     call reset_default_options(options)
     options%model = 1
     options%nlls_method = 1
     options%output_progress_vectors = .true.
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0) then
        write(*,*) 'Error: did not converge when output_progress_vectors = true'
        write(*,*) 'status = ', status%status
        no_errors_main = no_errors_main + 1
        status%status = 0
     end if
     options%output_progress_vectors = .false.

     ! Let's use a relative tr radius
     call reset_default_options(options)
     options%exact_second_derivatives = .true.
     options%relative_tr_radius = 1
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     options%relative_tr_radius = 0

     ! Let's use a weighted least squares method
     call reset_default_options(options)
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
     ! and the same with exact second derivatives (model = 2..4 )
     options%exact_second_derivatives = .true.
     do model = 2,4
        options%model = model
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
     end do
     deallocate(w)

     ! Let's do one run with non-exact second derivatives 
     call reset_default_options(options)
     options%exact_second_derivatives = .false.
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if

     ! And the same with model = 2
     call reset_default_options(options)
     options%model = 2
     options%exact_second_derivatives = .false.
     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        no_errors_main = no_errors_main + 1
     end if
     
     ! Let's get the method to switch to second-order
     call reset_default_options(options)
     options%exact_second_derivatives = .true.
     options%stop_g_absolute = 1e-14
     options%stop_g_relative = 1e-14
     options%print_level = 3
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
     call reset_default_options(options)
     options%model = 2
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
     call reset_default_options(options)
     options%nlls_method = 3125
     call solve_basic(X,params,options,status)
     if ( status%status .ne. ERROR%UNSUPPORTED_METHOD ) then 
        write(*,*) 'Error: unsupported method passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! test for unsupported tr strategy
     call reset_default_options(options)
     options%tr_update_strategy = 323
     call solve_basic(X,params,options,status)
     if ( status%status .ne. ERROR%BAD_TR_STRATEGY ) then 
        write(*,*) 'Error: unsupported TR strategy passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     
     ! three tests for incorrect returns from eval_f/J/H
     call reset_default_options(options)
     n = 2
     m = 67
     options%exact_second_derivatives = .true.
     do i = 1,3       
        X = [1.0, 2.0]
        select case (i)
        case (1)
           call nlls_solve(n, m, X,                         &
                eval_F_error, eval_J, eval_H, params,  &
                options, status )   
        case (2)
           call nlls_solve(n, m, X,                         &
                eval_F, eval_J_error, eval_H, params,  &
                options, status )   
        case (3)
           call nlls_solve(n, m, X,                         &
                eval_F, eval_J, eval_H_error, params,  &
                options, status )   
        end select
     end do
     if ( status%status .ne. ERROR%EVALUATION ) then 
        write(*,*) 'Error: error return from eval_x not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     if (no_errors_main == 0) then
        write(*,*) '*** All (main) tests passed successfully! ***'
     else
        write(*,*) 'There were ', no_errors_main,' errors'
     end if

     deallocate(model_to_test)
     deallocate(x)
     deallocate(params%x_values, params%y_values)
     
  end if

 
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
     

     call dogleg_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call generate_scaling_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call aint_tr_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call more_sorensen_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call evaluate_model_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call solve_galahad_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call solve_newton_tensor_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call all_eig_symm_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call solve_LLS_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call findbeta_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call calculate_rho_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call update_trust_region_radius_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call test_convergence_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call mult_J_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call mult_Jt_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails
     
     call switch_to_quasi_newton_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call solve_spd_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call solve_general_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call matmult_inner_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call matmult_outer_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call outer_product_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call min_eig_symm_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call max_eig_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call shift_matrix_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call get_svd_J_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails   
     
     call error_message_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails


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
