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
  real(wp), allocatable :: w(:),x(:),y(:),z(:)
  real(wp), allocatable :: A(:,:), B(:,:), C(:,:)
  real(wp), allocatable :: results(:)
  real(wp) :: alpha, beta, gamma, delta
  integer :: m, n, i, no_errors_helpers, no_errors_main, info
  integer :: total_errors
  integer :: nlls_method, model, tr_update
  logical :: test_all, test_subs

  integer :: number_of_models
  integer, allocatable :: model_to_test(:)

  type( NLLS_workspace ) :: work

  open(unit = 17, file="nlls_test.out")
  options%error = 17
  options%out   = 17 

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
          
     number_of_models = 7
     allocate(model_to_test(number_of_models))
     model_to_test = (/ 0, 1, 2, 3, 7, 8, 9 /)

     
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

              options%print_level = 0
              options%nlls_method = nlls_method
              options%tr_update_strategy = tr_update
              options%model = model_to_test(model)
              
              call nlls_solve(n, m, X,                         &
                   eval_F, eval_J, eval_H, params,  &
                   options, status )
              if (( nlls_method == 1).and.( options%model > 1)) then
                 if ( status%status .ne. -3 ) then
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
                 no_errors_main = no_errors_main + 1
              end if

           end do
        end do
     end do

     ! now let's check errors on the parameters passed to the routine...
     
     ! test n > m
     n = 100
     m = 3
     call  nlls_solve(n, m, X,                         &
                    eval_F, eval_J, eval_H, params,  &
                    options, status )
     if (status%status .ne. -800) then
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
        write(*,*) 'Error: unsupported TR strategy passedd and not caught'
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
     

     !! calculate step...
     ! ** TODO **

     !! dogleg 
     ! ** TODO **

     !! aint_tr
     ! ** TODO ** 

     !! more_sorensen
     ! ** TODO ** 

     !! solve_dtrs
     ! ** TODO ** 
     options%nlls_method = 4
     n = 2
     m = 5
     call setup_workspaces(work,n,m,options,info) 
     allocate(x(m*n),y(n),z(n**2),w(n))
     ! x -> J, y-> f, x -> hf, w-> d
     x = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp /)
     y = (/ 1.2_wp, 3.1_wp /)
     z = 1.0_wp

     alpha = 0.02_wp
     
     call solve_dtrs(x,y,z,n,m,alpha,w,& 
          work%calculate_step_ws%solve_dtrs_ws, &
          options,status)
     
     if ( status%status .ne. 0 ) then
        write(*,*) 'DTRS test failed, status = ', status%status
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     if ( abs(dot_product(w,w) - alpha**2) > 1e-12 ) then
        write(*,*) 'dtrs failed'
        write(*,*) 'Delta = ', alpha, '||d|| = ', dot_product(w,w)
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     deallocate(x,y,z,w)
     call remove_workspaces(work,options)

     !! solve_LLS 
     options%nlls_method = 1 ! dogleg
     call setup_workspaces(work,n,m,options,info) 
     

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
     call remove_workspaces(work, options)
     
     n = 100 
     m = 20
     allocate(x(n*m), y(m), z(n))     
     call setup_workspaces(work,n,m,options,info) 

     x = 1.0_wp
     x = 1.0_wp
     call solve_LLS(x,z,n,m,y,status, & 
          work%calculate_step_ws%dogleg_ws%solve_LLS_ws)
     if ( status%status .ne. ERROR%FROM_EXTERNAL ) then 
        write(*,*) 'solve_LLS test failed: wrong error message returned'
        write(*,*) 'status = ', status%status
        no_errors_helpers = no_errors_helpers+1
     end if
     
     deallocate(x,y,z)
     call remove_workspaces(work, options)
     options%nlls_method = 9 ! back to hybrid
     

     
     ! ** TODO **
     
         

     !------------!
     !! findbeta !!
     !------------!
     n = 3
     allocate(x(n),y(n),z(n))

     x = (/ 1.0, 2.0, 3.0 /) 
     y = (/ 2.0, 1.0, 1.0 /)

     call findbeta(x,y,1.0_wp,10.0_wp,alpha,info)

     if (info .ne. 0) then
        write(*,*) 'error -- findbeta did not work: info /= 0'
        no_errors_helpers = no_errors_helpers + 1
     else if ( ( norm2( x + alpha * y ) - 10.0_wp ) > 1e-12 ) then
        write(*,*) 'error -- findbeta did not work'
        write(*,*) '|| x + beta y|| = ', norm2( (x + alpha * y)-10.0_wp)
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
     call calculate_rho(alpha, beta, gamma, delta)
     if ( abs(delta - 3.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 3.0, got ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     
     ! now, let's check one is returned if alpha = beta
     beta = 2.0_wp
     call calculate_rho(alpha,beta,gamma, delta)
     if (abs(delta - 1.0_wp) > 1e-10) then
        write(*,*) 'Unexpected answer from calculate_rho'
        write(*,*) 'Expected 1.0, got ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     beta = 1.0_wp

     ! finally, check that 1 is returned if denominator = 0
     gamma = 2.0_wp
     call calculate_rho(alpha,beta,gamma, delta)
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
     
     delta = 100.0_wp ! Delta
     beta = 2.0_wp ! nu
     i = 3 ! p
     ! alpha = rho

     options%tr_update_strategy = 1
     ! let's go through the options
     
     options%eta_success_but_reduce = 0.25_wp
     options%eta_very_successful = 0.75_wp
     options%eta_too_successful = 2.0_wp

     ! check if rho reduced...
     alpha = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( delta >= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp
     
     ! check if rho stays the same...
     alpha = (options%eta_success_but_reduce + options%eta_very_successful) / 2
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( abs(delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp

     ! check if rho increases...
     alpha = (options%eta_very_successful + options%eta_too_successful) / 2
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( delta <= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp

     
     ! check if rho stays the same because too successful...
     alpha = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( abs(delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp

     ! now check for NaNs...HOW to do this in a non-compiler dependent way!?!?

     !! now, let's check the other option....
     options%tr_update_strategy = 1
     
     ! check if rho increases...
     alpha = (options%eta_very_successful + options%eta_too_successful) / 2
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( delta <= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not incease: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp
     
     ! check if rho stays the same because too successful...
     alpha = options%eta_too_successful + 1.0_wp
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( abs(delta - 100_wp) > 1e-12 ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not stay the same: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp

     alpha = options%eta_success_but_reduce - 0.5_wp
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( delta >= 100_wp ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Delta did not decrease as expected: delta = ', delta
        no_errors_helpers = no_errors_helpers + 1
     end if
     delta = 100_wp
     
     ! again...NaN test should go here!!!

     !Finally, check the error cases...
     
     options%tr_update_strategy = 18
     call update_trust_region_radius(alpha,options,status,delta,beta,i)
     if ( status%status .ne. ERROR%BAD_TR_STRATEGY ) then
        write(*,*) 'Unexpected answer from update_trust_region_radius'
        write(*,*) 'Error returned is = ', status%status, ', expected ',ERROR%BAD_TR_STRATEGY
        no_errors_helpers = no_errors_helpers + 1
     end if
     status%status = 0
     delta = 100_wp

     !! test_convergence
     ! todo
     
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
     call setup_workspaces(work,2,2,options,info) 
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

     call setup_workspaces(work,2,2,options,info) 
     allocate(A(n,n),x(2),y(2),z(2))

     A = reshape((/ 4.0, 1.0, 2.0, 2.0 /),shape(A))
     z = (/ 1.0, 1.0 /)
     y = (/ 6.0, 3.0 /)

     call solve_general(A,y,x,n,status,& 
          work%calculate_step_ws%AINT_tr_ws%solve_general_ws)
     if (status%status .ne. 0) then
        write(*,*) 'Error: info = ', info, ' returned from solve_general'
        no_errors_helpers = no_errors_helpers + 1
     else if (norm2(x-z) > 1e-12) then
        write(*,*) 'Error: incorrect value returned from solve_general'
        no_errors_helpers = no_errors_helpers + 1
     end if

     deallocate(A,x,y,z)
     call remove_workspaces(work,options)
     
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
        call setup_workspaces(work,n,m,options,info) 
        
        options%nlls_method = 2 ! revert...
        
        A = reshape( (/-5.0,  1.0, 0.0, 0.0, &
          1.0, -5.0, 0.0, 0.0, &
          0.0,  0.0, 4.0, 2.0, & 
          0.0,  0.0, 2.0, 4.0/), shape(A))

        call min_eig_symm(A,n,alpha,x,options,status, & 
             work%calculate_step_ws%more_sorensen_ws%min_eig_symm_ws)


        if ( (abs( alpha + 6.0 ) > 1e-12).or.(status%status .ne. 0) ) then
           write(*,*) 'error :: max_eig test failed -- wrong eig found'
           no_errors_helpers = no_errors_helpers + 1 
        elseif ( norm2(matmul(A,x) - alpha*x) > 1e-12 ) then
           write(*,*) 'error :: max_eig test failed -- not an eigenvector'
           no_errors_helpers = no_errors_helpers + 1
        end if
        
        call remove_workspaces(work,options)
        
     end do
          
     deallocate(A,x)

     !-----------!
     !! max_eig !!
     !-----------!
     n = 4
     m = 4
     ! make sure max_eig gets called

     allocate(x(n),A(n,n), B(n,n))
     call setup_workspaces(work,n,m,options,info) 
     
     A = reshape( (/1.0, 2.0, 3.0, 4.0, &
          2.0, 4.0, 6.0, 8.0, &
          3.0, 6.0, 9.0, 12.0, & 
          4.0, 8.0, 12.0, 16.0/), shape(A))
     B = 0.0_wp
     do i = 1,n
        B(i,i) = real(i,wp)
     end do
     
     
     call max_eig(A,B,n,alpha,x,C,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)

     if ( (abs( alpha - 10.0 ) > 1e-12).or.(status%status .ne. 0) ) then
        write(*,*) 'error :: max_eig test failed'
        no_errors_helpers = no_errors_helpers + 1 
     end if

     deallocate(A,B,x)
     call remove_workspaces(work,options)

     ! check the 'hard' case...
     n = 4
     m = 4
     allocate(x(n),A(n,n), B(n,n))
     call setup_workspaces(work,2,2,options,info) 

     A = 0.0_wp  
     A(3,1) = 1.0_wp; A(4,1) = 2.0_wp; A(3,2) = 3.0_wp; A(4,2) = 4.0_wp
     A(1,3) = A(3,1); A(1,4) = A(4,1); A(2,3) = A(3,2); A(2,4) = A(4,2)
     B = A
     A(1,1) = 1.0_wp; A(2,2) = 1.0_wp

     call max_eig(A,B,n,alpha,x,C,options,status, & 
                  work%calculate_step_ws%AINT_tr_ws%max_eig_ws)

     if (.not. allocated(C)) then ! check C returned 
        write(*,*) 'error :: hard case of max_eig test failed - C not returned'
        no_errors_helpers = no_errors_helpers + 1 
     else
        allocate(y(2))
        y = shape(C)
        if ((y(1) .ne. 2) .or. (y(2) .ne. n)) then
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

     deallocate(A,B,C,x,y)
     if (allocated(results)) deallocate(results)
     call remove_workspaces(work,options)

     
     
     call setup_workspaces(work,1,1,options,info)  !todo: deallocation routine
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
     call remove_workspaces(work,options)

     !! shift_matrix 
     ! todo
     
     !! get_svd_J 
     ! Todo

     ! Report back results....

     if (no_errors_helpers == 0) then
        write(*,*) '*** All (helper) tests passed successfully! ***'
     else
        write(*,*) 'There were ', no_errors_helpers,' errors'
     end if

  end if

  
close(unit = 17)
!no_errors_helpers = 1
 if (no_errors_helpers + no_errors_main == 0) then
    write(*,*) ' '
    write(*,*) '**************************************'
    write(*,*) '*** All tests passed successfully! ***'
    write(*,*) '**************************************'
!    stop 0    ! needed for talking with ctest
 else 
!    stop 1    ! needed for talking with ctest
  end if
  


end program nlls_test
