! Copyright (c) 2019, The Numerical Algorithms Group Ltd (NAG)
! All rights reserved.
! Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

#include "../src/preprocessor.FPP"

program nlls_test

  use MODULE_PREC(unit_test_mod)

  implicit none

  type( NLLS_inform )  :: status
  type( NLLS_options ) :: options
  type( user_type ), target :: params
!  type( user_box_type ), target :: params_box
  real(wp), allocatable :: w(:),x(:),blx(:),bux(:)
  real(wp), allocatable :: resvec(:)
  real(wp) :: resvec_error, tol
  integer :: m, n, i, no_errors_helpers, no_errors_main
  integer :: nlls_method, model, tr_update, inner_method
  logical :: test_all, test_subs, oki
  integer :: fails, exp_status

  integer :: number_of_models
  integer, allocatable :: model_to_test(:)
!  character*40 :: details

  options%out   = 17
  options%print_level = 0
  open(unit = options%out, file="nlls_test.out")

  test_all = .true.
  test_subs = .true.
  exp_status = 0
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
     options%print_level = 5

     do tr_update = 1,2
        do nlls_method = 1,4
           do model = 1,number_of_models

              call reset_default_options(options)
              options%print_options = .True.
              options%nlls_method = nlls_method
              options%tr_update_strategy = tr_update
              options%model = model_to_test(model)
              options%exact_second_derivatives = .true.
              options%output_progress_vectors = .true.
              call print_line(options%out)
              write(options%out,*) "tr_update_strategy = ", options%tr_update_strategy
              write(options%out,*) "nlls_method        = ", options%nlls_method
              write(options%out,*) "model              = ", options%model
              call print_line(options%out)
              if (nlls_method == 4) then
                 do inner_method = 1,3
                    ! check the tests with c and fortran jacobians
                    ! pass individually, and give consistent results.
                    options%inner_method = inner_method
                    call c_fortran_tests(options,no_errors_main)
                 end do
              else
                 call c_fortran_tests(options,no_errors_main)
              end if
           end do
        end do
     end do


     ! dogleg, no fallback
     call reset_default_options(options)
     options%nlls_method = 1 ! dogleg
     options%model = 2 ! newton
     options%allow_fallback_method = .false.
     options%print_level = 1

     call print_line(options%out)
     write(options%out,*) "dogleg, model = 2, no fallback"
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_DOGLEG_MODEL ) then
        write(*,*) 'incorrect error return from nlls_solve:'
        write(*,*) 'NLLS_METHOD = ', nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'status returned = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! now, let's test the regularization method
     call reset_default_options(options)
     options%type_of_method = 2 ! regularization
     options%inner_method = 2
     options%print_level = 1
     options%exact_second_derivatives = .true.

     call print_line(options%out)
     write(options%out,*) "type_of_method = ", options%type_of_method
     call print_line(options%out)


     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: regularization'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! now, let's run through all the print options...
     do i = 1,6

        call reset_default_options(options)
        options%print_level = i

        call print_line(options%out)
        write(options%out,*) "print_level = ", options%print_level
        call print_line(options%out)


        call solve_basic(X,params,options,status)
        if ( i == 6) then
           if ( status%status .ne. NLLS_ERROR_PRINT_LEVEL) then
              write(*,*) 'print test: expected error not returned'
              write(*,*) 'info%status = ', status%status
              no_errors_main = no_errors_main + 1
           end if
        elseif (status%status .ne. 0) then
           write(*,*) 'nlls_solve failed to converge: print test'
           write(*,*) 'info%status = ', status%status
           no_errors_main = no_errors_main + 1
        end if
     end do

     ! and the print options with regularization....
     call reset_default_options(options)
     options%type_of_method = 2 ! regularization
     options%inner_method = 2
     options%exact_second_derivatives = .true.

     do i = 1,6

        options%print_level = i

        call print_line(options%out)
        write(options%out,*) "type_of_method = ", options%type_of_method
        write(options%out,*) "print_level = ", options%print_level
        call print_line(options%out)
        call solve_basic(X,params,options,status)
        if ( i == 6) then
           if ( status%status .ne. NLLS_ERROR_PRINT_LEVEL) then
              write(*,*) 'print test: expected error not returned'
              write(*,*) 'info%status = ', status%status
              no_errors_main = no_errors_main + 1
           end if
        elseif (status%status .ne. 0) then
           write(*,*) 'nlls_solve failed to converge: print test'
           write(*,*) 'info%status = ', status%status
           no_errors_main = no_errors_main + 1
        end if

     end do

     options%print_level = 1


     do i = 1,2
        ! Let's do a test where the regularization weight is non-zero
        call reset_default_options(options)
        options%regularization = i
        options%regularization_term = 1e-2
        options%regularization_power = 2.0_wp

        call print_line(options%out)
        write(options%out,*) "Regularization power is two, weight is non-zero"
        call print_line(options%out)

        call solve_basic(X,params,options,status)
        if ( status%status .ne. 0 ) then
           write(*,*) 'nlls_solve failed to converge: non-zero regularization weight'
           write(*,*) 'NLLS_METHOD = ', options%nlls_method
           write(*,*) 'MODEL = ', options%model
           write(*,*) 'info%status = ', status%status
           no_errors_main = no_errors_main + 1
        end if

        call print_line(options%out)
        write(options%out,*) "Regularization power is three, weight is non-zero"
        call print_line(options%out)

        ! and, the same test with a regularization power of three:
        call reset_default_options(options)
        options%regularization = i
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

        call print_line(options%out)
        write(options%out,*) "Regularization, model = 2"
        call print_line(options%out)

        ! now let's get regularization with model = 2
        call reset_default_options(options)
        options%regularization = i
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
     end do

     ! now test optimal reg_order
     call reset_default_options(options)
     options%type_of_method = 2
     options%model = 1
     options%reg_order = -1.0

     call print_line(options%out)
     write(options%out,*) "Optimal regularization order"
     call print_line(options%out)


     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: negative reg_order'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! now test optimal reg_order, hybrid model
     call reset_default_options(options)
     options%type_of_method = 2
     options%model = 3
     options%reg_order = -1.0


     call print_line(options%out)
     write(options%out,*) "Optimal regularization order with hybrid model"
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: negative reg_order'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! now test home rolled regularization
     call reset_default_options(options)
     options%type_of_method = 2  ! regularization
     options%model = 1 ! Gauss-Newton
     options%nlls_method = 3
     options%reg_order = 2.0
     options%print_level = 5

     call print_line(options%out)
     write(options%out,*) "Regularization"
     write(options%out,*) "model = ", options%model
     write(options%out,*) "nlls_method = ", options%nlls_method
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: negative reg_order'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! home-rolled MS, no eigenvalues
     call reset_default_options(options)
     options%model = 1 ! Gauss-Newton
     options%nlls_method = 3
     options%use_ews_subproblem = .false.
     options%print_level = 5

     call print_line(options%out)
     write(options%out,*) "Trust region"
     write(options%out,*) "model = ", options%model
     write(options%out,*) "nlls_method = ", options%nlls_method
     write(options%out,*) "use_ews_subproblem = ", options%use_ews_subproblem
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. 0 ) then
        write(*,*) 'no eigenvalue MS failed to converge'
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
     options%print_level = 4

     call print_line(options%out)
     write(options%out,*) "Regularization"
     write(options%out,*) "model = ", options%model
     call print_line(options%out)

     do inner_method = 1,4
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

     ! let's try one with HP included...
     call print_line(options%out)
     write(options%out,*) "Pass in eval_HP"
     call print_line(options%out)

     options%inner_method = 2
     n = 2
     m = 67
     X = [1.0, 2.0]
     call nlls_solve(n, m, X,                        &
          eval_F,eval_J,eval_H, params,   &
          options, status, eval_HP=eval_HP)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: tensor model, eval_HP'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if


     ! and HP + weights....
     options%inner_method = 2
     n = 2
     m = 67
     X = [1.0, 2.0]
     allocate(w(m))
     do i = 1, m
        w(i) = 2.0
     end do

     call print_line(options%out)
     write(options%out,*) "Pass in eval_HP and weights"
     call print_line(options%out)


     call nlls_solve(n, m, X,                        &
          eval_F,eval_J,eval_H, params,   &
          options, status, weights=w, eval_HP=eval_HP)
     if ( status%status .ne. 0 ) then
        write(*,*) 'nlls_solve failed to converge: tensor model, eval_HP, weights'
        write(*,*) 'NLLS_METHOD = ', options%nlls_method
        write(*,*) 'MODEL = ', options%model
        write(*,*) 'info%status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     deallocate(w)

     ! now, let's get an error return...
     call reset_default_options(options)
     options%model = 4
     options%exact_second_derivatives = .false.
     call print_line(options%out)
     write(options%out,*) "model = ", options%model
     write(options%out,*) "exact_second_derivatives = ", options%exact_second_derivatives
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_NO_SECOND_DERIVATIVES ) then
        write(*,*) 'expected error return', NLLS_ERROR_NO_SECOND_DERIVATIVES,' but'
        write(*,*) 'got ', status%status
        no_errors_main = no_errors_main + 1
     end if

     ! Let's get to maxits
     call reset_default_options(options)
     options%type_of_method = 1
     options%maxit = 5
     options%model = 1
     options%nlls_method = 1

     call print_line(options%out)
     write(options%out,*) "Reach maxits"
     call print_line(options%out)

     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_MAXITS) then
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

     call print_line(options%out)
     write(options%out,*) "output_progress_vectors = ", options%output_progress_vectors
     call print_line(options%out)

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
     call print_line(options%out)
     write(options%out,*) "relative_tr_radius = ", options%relative_tr_radius
     call print_line(options%out)

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

     call print_line(options%out)
     write(options%out,*) "Pass in weights"
     call print_line(options%out)

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

        call print_line(options%out)
        write(options%out,*) "Pass in weights"
        write(options%out,*) "model = ", options%model
        call print_line(options%out)
        call nlls_solve(n, m, X,                         &
             eval_F, eval_J, eval_H, params,  &
             options, status, weights=w )
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

     call print_line(options%out)
     write(options%out,*) "exact_second_derivatives = ", options%exact_second_derivatives
     call print_line(options%out)


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
     call print_line(options%out)
     write(options%out,*) "exact_second_derivatives = ", options%exact_second_derivatives
     write(options%out,*) "model = ", options%model
     call print_line(options%out)

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
     call print_line(options%out)
     write(options%out,*) "switch to second order"
     call print_line(options%out)

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
     call print_line(options%out)
     write(options%out,*) "||f|| = 0"
     call print_line(options%out)

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

     ! test for c-based Jabobians
     ! run fortran based and c based, and check the resvecs are the same
     call reset_default_options(options)
     ! first, let's get a standard output...
     call print_line(options%out)
     write(options%out,*) "C-based Jacobians: get standard output"
     call print_line(options%out)

     n = 2
     m = 67
     X =[1.0, 2.0]
     options%output_progress_vectors = .true.

     call nlls_solve(n, m, X,                         &
                     eval_F, eval_J, eval_H, params,  &
                     options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'solve failed'
        no_errors_main = no_errors_main + 1
     end if
     ! save the resvec, so that we can compare later
     allocate(resvec(status%iter))
     resvec(1:status%iter) = status%resvec(1:status%iter)

     ! now run with a row-major ordered Jacobian
     X = [1.0, 2.0]
     options%Fortran_Jacobian = .false.
     call print_line(options%out)
     write(options%out,*) "C-based Jacobians: Fortran_Jacobian = ", options%Fortran_Jacobian
     call print_line(options%out)

     call nlls_solve(n, m, X,                         &
                     eval_F, eval_J_c, eval_H, params,  &
                     options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'solve failed'
        no_errors_main = no_errors_main + 1
     end if
     resvec(:) = resvec(:) - status%resvec(1:size(resvec))
     resvec_error = dot_product( resvec(:),resvec(:) )
     if ( resvec_error > 1e-14 ) then
        write(*,*) 'C Jacobian test failed: resvec error = ', resvec_error
        write(*,*) 'resvec = ', resvec
        no_errors_main = no_errors_main + 1
     end if


     deallocate(resvec)

     ! now run with a row-major ordered Jacobian
     ! and also with a relative tr radius
     X = [1.0, 2.0]
     options%Fortran_Jacobian = .false.
     options%relative_tr_radius = 1
     call print_line(options%out)
     write(options%out,*) "C-based Jacobians: Fortran_Jacobian = ", options%Fortran_Jacobian
     write(options%out,*) "relative_tr_radius = ", options%relative_tr_radius
     call print_line(options%out)

     call nlls_solve(n, m, X,                         &
                     eval_F, eval_J_c, eval_H, params,  &
                     options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'C Jacobian, relative_tr_radius test failed'
        no_errors_main = no_errors_main + 1
     end if

     ! three tests for incorrect returns from eval_f/J/H
     call reset_default_options(options)
     n = 2
     m = 67
     options%exact_second_derivatives = .true.
     do i = 1,5
        X = [1.0, 2.0]
        select case (i)
        case (1)
           call print_line(options%out)
           write(options%out,*) "Error in eval_F"
           call print_line(options%out)
           exp_status = NLLS_ERROR_INITIAL_GUESS

           call nlls_solve(n, m, X,                         &
                eval_F_error, eval_J, eval_H, params,  &
                options, status )
        case (2)
           call print_line(options%out)
           write(options%out,*) "Error in eval_J"
           call print_line(options%out)
           exp_status = NLLS_ERROR_INITIAL_GUESS

           call nlls_solve(n, m, X,                         &
                eval_F, eval_J_error, eval_H, params,  &
                options, status )
        case (3)
           call print_line(options%out)
           write(options%out,*) "Error in eval_HF"
           call print_line(options%out)
           exp_status = NLLS_ERROR_EVALUATION

           call nlls_solve(n, m, X,                         &
                eval_F, eval_J, eval_H_error, params,  &
                options, status )
        case (4)
           options%model = 2
           call print_line(options%out)
           write(options%out,*) "Error in eval_HF"
           write(options%out,*) "model = ", options%model
           call print_line(options%out)
           exp_status = NLLS_ERROR_EVALUATION

           call nlls_solve(n, m, X,                         &
                eval_F, eval_J, eval_H_error, params,  &
                options, status )
        case (5)
           call print_line(options%out)
           write(options%out,*) "Error in eval_HF"
           write(options%out,*) "model = ", options%model
           call print_line(options%out)
           exp_status = NLLS_ERROR_INITIAL_GUESS

           options%model = 4
           call nlls_solve(n, m, X,                         &
                eval_F, eval_J, eval_H_error, params,  &
                options, status )
        end select
        if ( status%status .ne. exp_status ) then
           write(*,*) 'Error: error return from eval_x not caught'
           no_errors_main = no_errors_main + 1
        end if
     end do
     status%status = 0


     ! three tests for incorrect returns from eval_f/J/H
     ! after the first case
     call reset_default_options(options)
     n = 2
     m = 67
!     options%exact_second_derivatives = .true.
     do i = 1,3
        X = [1.0, 2.0]
        params%iter = 0
        select case (i)
        case (1)
           call print_line(options%out)
           write(options%out,*) "Error in eval_F at iteration 2"
           call print_line(options%out)
           call nlls_solve(n, m, X,                         &
                eval_F_one_error, eval_J, eval_H, params,  &
                options, status )
        case (2)
           call print_line(options%out)
           write(options%out,*) "Error in eval_J at iteration 2"
           call print_line(options%out)
           call nlls_solve(n, m, X,                        &
                eval_F, eval_J_one_error, eval_H, params,  &
                options, status )
        case (3)
           call print_line(options%out)
           write(options%out,*) "Error in eval_HF at iteration 2"
           call print_line(options%out)
           call nlls_solve(n, m, X,                        &
                eval_F, eval_J, eval_H_one_error, params,  &
                options, status )
        end select
        if ( status%status .ne. 0 ) then
           write(*,*) 'Error: single error return from eval_x should have worked'
           write(*,*) '       but status = ', status%status, ' returned'
           no_errors_main = no_errors_main + 1
        end if
     end do
     status%status = 0


     ! tests for too many reductions of tr
     call reset_default_options(options)
     n = 2
     m = 67
!     options%exact_second_derivatives = .true.
     X = [1.0, 2.0]
     params%iter = 0
     call print_line(options%out)
     write(options%out,*) "Too many TR radius reductions"
     call print_line(options%out)
     call nlls_solve(n, m, X,                         &
          eval_F_allbutone_error, eval_J, eval_H, params,  &
          options, status )
     if ( status%status .ne. NLLS_ERROR_MAX_TR_REDUCTIONS ) then
        write(*,*) 'Error: expected to many reductions error'
        write(*,*) '       but status = ', status%status, ' returned'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! tests for too large f(x)
     call reset_default_options(options)
     n = 2
     m = 67
!     options%exact_second_derivatives = .true.
     X = [1.0, 2.0]
     params%iter = 0
     options%print_level = 5
     call print_line(options%out)
     write(options%out,*) "Large JtF, defaults"
     call print_line(options%out)
     call nlls_solve(n, m, X,                         &
          eval_F_large, eval_J_large, eval_H, params,  &
          options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'Error: Large JtF test failed, but pass expected'
        write(*,*) '       status = ', status%status, ' returned'
        write(*,*) status%error_message
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! tests for too large f(x) with exact second derivatives
     call reset_default_options(options)
     n = 2
     m = 67
     options%exact_second_derivatives = .true.
     X = [1.0, 2.0]
     params%iter = 0
     options%print_level = 3
     options%print_options = .true.

     call print_line(options%out)
     write(options%out,*) "Large JtF (exact H)"
     call print_line(options%out)
     call nlls_solve(n, m, X,                         &
          eval_F, eval_J, eval_H, params,  &
          options, status )

     X = [1.0, 2.0]
     call nlls_solve(n, m, X,                         &
          eval_F_large, eval_J_large, eval_H, params,  &
          options, status )
     if ( status%status .ne. 0 ) then
        write(*,*) 'Error: Large JtF test (exact H) failed, but pass expected'
        write(*,*) '       status = ', status%status, ' returned'
        write(*,*) status%error_message
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0


     ! tests for too large f(x) with weights
     call reset_default_options(options)
     n = 2
     m = 67
     X = [1.0, 2.0]

     allocate(w(m))
     do i = 1, m
        w(i) = 2.0
     end do

     params%iter = 0
     options%print_level = 2
     call print_line(options%out)
     write(options%out,*) "Large JtF (with weights)"
     call print_line(options%out)
     call nlls_solve(n, m, X,                         &
          eval_F_large, eval_J_large, eval_H, params,  &
          options, status, weights=w )
     if ( status%status .ne. 0 ) then
        write(*,*) 'Error: Large JtF (with weights) test failed, but pass expected'
        write(*,*) '       status = ', status%status, ' returned'
        write(*,*) status%error_message
        no_errors_main = no_errors_main + 1
     end if
     deallocate(w)
     status%status = 0



     ! now let's check errors on the parameters passed to the routine...

     call print_line(options%out)
     write(options%out,*) "Parameter errors"
     call print_line(options%out)

     options%print_level = 3

    ! test for unsupported method
     call reset_default_options(options)
     options%nlls_method = 3125
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_UNSUPPORTED_METHOD ) then
        write(*,*) 'Error: unsupported method passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

    ! test for unsupported inner_method
     call reset_default_options(options)
     options%model = 4
     options%inner_method = 3125
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_WRONG_INNER_METHOD ) then
        write(*,*) 'Error: wrong inner method passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! test for unsupported method, with type_of_method = 2
     call reset_default_options(options)
     options%type_of_method = 2
     options%nlls_method = 3125
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_UNSUPPORTED_METHOD ) then
        write(*,*) 'Error: unsupported method (nlls_method=2) passed and not caught'
        write(*,*) 'status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! test for unsupported type_of_method
     call reset_default_options(options)
     options%type_of_method = 2343
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_UNSUPPORTED_TYPE_METHOD ) then
        write(*,*) 'Error: unsupported type_of_method passed and not caught'
        write(*,*) 'status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! test for unsupported tr strategy
     call reset_default_options(options)
     options%tr_update_strategy = 323
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_BAD_TR_STRATEGY ) then
        write(*,*) 'Error: unsupported TR strategy passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! Test for initial point being a solution
     ! Two cases:
     ! 1. x is solution point
     ! 2. x is active and projected gradient is zero
     deallocate(params%x_values, params%y_values)
     call generate_data_example_box(params)
     call reset_default_options(options)
     options%maxit = 2
     options%print_level=0
     if (wp == lp) then
        ! relax convergence criteria
        options%stop_g_absolute = 1.0e-4
     end if
     x(:) = (/0.3199787042575630E+00, 0.2752509146444680E-01/)
     call solve_basic(X,params,options,status,.True.)
     if ( .Not. (status%status == 0 .And. status%iter == 0) ) then
        write(*,*) 'Error: x0 solution but not caught'
        no_errors_main = no_errors_main + 1
     end if

     call reset_default_options(options)
     options%maxit = 2
     options%print_level=0
     if (wp == lp) then
        ! relax convergence criteria
        options%stop_g_absolute = 1.0e-4
     end if
     x(:) = (/0.3199787042575630E+00, 0.2752509146444680E-01/)
     Allocate(blx(n), bux(n))
     blx(1:n) = -1.0
     bux(1:n) = blx(1:n)
     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( .Not. (status%status == 0 .And. status%iter == 0 .And.               &
       status%norm_g==0.0) ) then
        write(*,*) 'Error: Proj grd at x0 is zero, but not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! Unsupported Line Search
     call reset_default_options(options)
     options%box_linesearch_type = 3
     options%print_level = 5
     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( status%status .ne. NLLS_ERROR_UNSUPPORTED_LINESEARCH ) then
        write(*,*) 'Error: unsupported Linesearch type passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! Unsupported Print Level
     call reset_default_options(options)
     options%print_level = 7
     call solve_basic(X,params,options,status)
     if ( status%status .ne. NLLS_ERROR_PRINT_LEVEL ) then
        write(*,*) 'Error: unsupported print level passed and not caught'
        no_errors_main = no_errors_main + 1
     end if
     options%print_level = 0
     status%status = 0

     ! Bad box constraints
     call reset_default_options(options)
     blx(1:n) = 1.0
     bux(1:n) = -1.0
     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( status%status /= NLLS_ERROR_BAD_BOX_BOUNDS ) then
        write(*,*) 'Error: Illegal box, but not caught.  Status = ', status%status
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! Projected Gradient Linesearch recovery error after LS D-S failure
     call reset_default_options(options)
     x(1) = 1.3_wp
     x(2) = -2.0_wp
     blx(1:n) = (/1.2_wp,-10.0_wp/)
     bux(1:n) = (/1.2_wp,10.0_wp/)
     call nlls_solve(n, m, x,                         &
          eval_F_pg, eval_J, eval_H, params,  &
          options, status,                 &
          lower_bounds=blx, upper_bounds=bux )
     if ( status%status /= NLLS_ERROR_PG_STEP ) then
        write(*,*) 'Error: PG step failed, but not caught'
        write(*,*) 'status = ', status%status, ' returned'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0

     ! Line-search HZ supported but not yet available LINE SEARCH ERROR
     call reset_default_options(options)
     x(1) = 1.3_wp
     x(2) = -2.0_wp
     blx(1:n) = (/1.2_wp,-10.0_wp/)
     bux(1:n) = (/1.2_wp,10.0_wp/)
     options%box_linesearch_type = 2
     options%maxit = 10
     call solve_basic(X,params,options,status,warm_start=.True.,blx=blx,bux=bux)
     if ( status%status /= NLLS_ERROR_UNSUPPORTED_LINESEARCH ) then
        write(*,*) 'Error: HZLS unsupported, but not caught.'
        write(*,*) 'status = ', status%status, ' returned'
        no_errors_main = no_errors_main + 1
     end if
     status%status = 0
     options%print_level = 0
     options%out = 17

     ! Excercise Print Level logic, this is not a test
     call reset_default_options(options)
     blx(1:n) = (/1.2_wp,-10.0_wp/)
     bux(1:n) = (/1.2_wp,10.0_wp/)
     Do i = 1, 5
       options%print_level = i
       Write(options%out, '(80(''=''))')
       Write(options%out, *) 'Exercising Print Level = ', options%print_level
       Write(options%out, '(80(''=''))')
       Call reset_default_options(options)
       options%print_header = 4
       options%print_options = .True.
       x(1) = 1.3_wp
       x(2) = -2.0_wp
       call solve_basic(X,params,options,status,blx=blx,bux=bux)
       status%status = 0
     End Do
     status%status = 0

     ! Check invalid weights vector
     deallocate(params%x_values, params%y_values)
     call generate_data_example_box(params)
     call reset_default_options(options)
     options%maxit = 1
     options%print_level=1
     x(:) = (/0.3199787042575630E+00, 0.2752509146444680E-01/)
     status%status = 0
     Allocate(w(m))
     w(1:m) = -1.0
     call solve_basic(X,params,options,status,weights=w)
     if ( status%status /= NLLS_ERROR_BAD_WEIGHTS ) then
        write(*,*) 'Error: Solver did not detect illegal weights'
        no_errors_main = no_errors_main + 1
     end if

     ! ========================================================================
     ! Finite-Differences Unit-Tests
     ! Hit bound on either side
     deallocate(params%x_values, params%y_values)
     call generate_data_example_box(params)
     call reset_default_options(options)
     options%maxit = 2
     options%print_level = 2
     options%check_derivatives = 1
     if (wp == lp) then
        options%fd_step = 1.0e-3_wp
        options%derivative_test_tol = 5.0e-3_wp
     end if
     x(:) = (/0.0, 0.0/)
     status%status = 0
     blx(1:n) = 1.0
     bux(1:n) = 2.0

     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( .Not. ( status%status == -1 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if

     ! catch error in Jacobian:
     options%print_level = 1
     call nlls_solve(n, m, x,                         &
          eval_F, eval_J_bad, ral_nlls_eval_hf_dummy, params,  &
          options, status )
     if ( .Not. ( status%status == -19 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if
     options%print_level = 0

     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( .Not. ( status%status == -1 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if

     ! check again with weights
     w(1:m) = 1.5
     call solve_basic(X,params,options,status,blx=blx,bux=bux,weights=w)
     if ( .Not. ( status%status == -1 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if

     ! Wrong storage order
     options%Fortran_Jacobian = .False.
     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( .Not. ( status%status == -19) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if
     options%Fortran_Jacobian = .True.

     ! No search space -> can't check any: box too tight
     blx(1:n) = 1.0
     bux(1:n) = 1.0
     call solve_basic(X,params,options,status,blx=blx,bux=bux)
     if ( .Not. ( status%status == 0 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if

     ! Fail the derivative test (too stringent)
     options%derivative_test_tol = 1.0e-16_wp
     options%print_level = 1
     call solve_basic(X,params,options,status)
     if ( .Not. ( status%status == -19 ) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if
     options%derivative_test_tol = 1.0e-4_wp

     ! Solve the problem using FD (Fortran storage)
     options%Fortran_Jacobian = .True.
     options%print_level = 2
     options%check_derivatives = 0
     options%maxit = 100
     call solve_basic(X,params,options,status,use_fd=.True.)
     if (wp == np) then
        tol = 1.e-6_wp
     else
        ! Relax tolerance for lower precision
        tol = 5e-5_wp
     end if
     oki = abs(x(1) - 0.319978704257563_wp) <= tol .And. abs(x(2)-0.2752509146444680E-1_wp) < tol
     write(options%out,*) 'Solution (F): ', x(1:n)
     write(options%out,*) 'Expected: ', (/0.319978704257563_wp, 0.2752509146444680E-1_wp/)
     if ( .Not. ( status%status == 0 .And. oki) ) then
        write(*,*) 'Error: FD solve: unexpected status value or wrong solution'
        no_errors_main = no_errors_main + 1
     end if

     ! Solve the problem using FD (C storage)
     options%Fortran_Jacobian = .False.
     options%print_level = 3
     call solve_basic(X,params,options,status,use_fd=.True.)
     oki = abs(x(1) - 0.3199787) <= tol .And. abs(x(2)-0.02752509) < tol
     write(options%out,*) 'Solution (C): ', x(1:n)
     write(options%out,*) 'Expected: ', (/0.3199787042575630E+00, 0.2752509146444680E-01/)
     if ( .Not. ( status%status == 0 .And. oki) ) then
        write(*,*) 'Error: FD solve: unexpected status value or wrong solution'
        no_errors_main = no_errors_main + 1
     end if

     ! Solve with one fixed variable
     options%Fortran_Jacobian = .False.
     options%print_options = .True.
     blx(:) = (/0.0, 0.0/)
     bux(:) = (/0.0, 1.0/)
     if (wp == np) then
        tol = 1.e-6_wp
        ! options%box_gamma = 0.99999_wp
     else
        ! Relax tolerance for lower precision
        tol = 2e-3_wp
     end if
     call solve_basic(X,params,options,status,use_fd=.True.,blx=blx,bux=bux)
     oki = x(1) == 0.0 .And. abs(x(2)-0.923618046017) < tol
     write(options%out,*) 'Solution (C): ', x(1:n)
     write(options%out,*) 'Expected: ', (/0.0, 0.92361804601/)
     write(options%out,*) 'DIFF: ', abs(x(2)-0.92361804601), '(',tol,')'
     if ( .Not. ( status%status == 0 .And. oki) ) then
        write(*,*) 'Error: FD solve 3: unexpected status value or wrong solution'
        write(*,*) '       status = ', status%status, "(expected 0)"
        write(*,*) '       Solution (C): ', x(1:n)
        write(*,*) '       Expected: ', (/0.0, 0.92361804601/)
        write(*,*) '       DIFF: ', abs(x(2)-0.92361804601), '(',tol,')'
        no_errors_main = no_errors_main + 1
     end if

     Call reset_default_options(options)

     ! Solve with one active variable
     options%Fortran_Jacobian = .True.
     options%print_level = 3
     options%check_derivatives = 0
     blx(:) = (/0.32, 0.0/)
     bux(:) = (/1.0, 1.0/)
     if (wp == np) then
        tol = 1.e-6_wp
        options%maxit = 100
     else
        ! Relax tolerance for lower precision
        tol = 2.0e-5_wp
        options%fd_step = 9.0e-3_wp
        options%maxit = 200
     end if
     call solve_basic(X,params,options,status,use_fd=.True.,blx=blx,bux=bux)
     oki = abs(x(1) - 0.32) <= tol .And. abs(x(2)-2.7447762174312485E-2) < tol
     write(options%out,*) 'Solution (F): ', x(1:n)
     write(options%out,*) 'Expected: ', (/0.32, 2.7447762174312485E-2/)
     if ( .Not. ( status%status == 0 .And. oki) ) then
        write(*,*) 'Error: FD solve: unexpected status value or wrong solution'
        no_errors_main = no_errors_main + 1
     end if

     ! eval_F fails midway FD estimation - no recovery
     options%check_derivatives = 1
     params%iter = 0
     call nlls_solve(n, m, x,                         &
          eval_F_one_error, eval_J, ral_nlls_eval_hf_dummy, params,  &
          options, status )
     if ( .Not. ( status%status == -4 .And. status%external_return==2101) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if
     ! try now with rubbish eval_f residual - no recovery
     params%iter = 0
     call nlls_solve(n, m, x,                         &
          eval_F_one_NaN, eval_J, ral_nlls_eval_hf_dummy, params,  &
          options, status )
     if ( .Not. ( status%status == -4 .And. status%external_return==2031) ) then
        write(*,*) 'Error: check_derivatives: unexpected status value'
        no_errors_main = no_errors_main + 1
     end if
     options%check_derivatives = 0

     status%status = 0
     call reset_default_options(options)


     Write(options%out, '(80(''=''))')
     options%print_level = 0

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

     call trust_region_subproblem_tests(options,fails)
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

     call minus_solve_spd_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call minus_solve_general_tests(options,fails)
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

     call error_message_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call covariance_matrix_tests(options,fails)
     no_errors_helpers = no_errors_helpers + fails

     call evaltensor_J_tests(options,fails)
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
    write(*,*) ' '
    write(*,*) ' Error summary for all tests:'
    write(*,*) '  TOTAL:           ', no_errors_main + no_errors_helpers
    write(*,*) '    Main tests:    ', no_errors_main
    write(*,*) '    Helpers tests: ', no_errors_helpers
    write(*,*) ' '
    write(*,*) ' ** UNIT TESTS FAILED **'
    stop 1    ! needed for talking with ctest
  end if



end program nlls_test
