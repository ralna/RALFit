! ral_nlls_internal :: internal subroutines for ral_nlls

module ral_nlls_internal

  use RAL_NLLS_DTRS_double
  use RAL_NLLS_DRQS_double
  use ral_nlls_workspaces
  
  implicit none

  private

  integer, parameter :: wp = kind(1.0d0)
  integer, parameter :: long = selected_int_kind(8)
  real (kind = wp), parameter :: tenm3 = 1.0e-3
  real (kind = wp), parameter :: tenm5 = 1.0e-5
  real (kind = wp), parameter :: tenm8 = 1.0e-8
  real (kind = wp), parameter :: epsmch = epsilon(1.0_wp)
  real (kind = wp), parameter :: hundred = 100.0
  real (kind = wp), parameter :: ten = 10.0
  real (kind = wp), parameter :: point9 = 0.9
  real (kind = wp), parameter :: zero = 0.0
  real (kind = wp), parameter :: one = 1.0
  real (kind = wp), parameter :: two = 2.0
  real (kind = wp), parameter :: three = 3.0
  real (kind = wp), parameter :: half = 0.5
  real (kind = wp), parameter :: sixteenth = 0.0625
  
  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: f
       class(params_base_type), intent(in) :: params
     end subroutine eval_f_type
  end interface

  abstract interface
     subroutine eval_j_type(status, n, m, x, J, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: J
       class(params_base_type), intent(in) :: params
     end subroutine eval_j_type
  end interface

  abstract interface
     subroutine eval_hf_type(status, n, m, x, f, h, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m 
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: f
       double precision, dimension(*), intent(out) :: h
       class(params_base_type), intent(in) :: params
     end subroutine eval_hf_type
  end interface

    public :: nlls_solve, nlls_iterate, nlls_finalize, nlls_strerror
    public :: solve_galahad, findbeta, mult_j
    public :: mult_jt, solve_spd, solve_general, matmult_inner
    public :: matmult_outer, outer_product, min_eig_symm, max_eig, all_eig_symm
    public :: remove_workspaces, setup_workspaces
    public :: get_svd_j, calculate_step, evaluate_model
    public :: update_trust_region_radius, apply_second_order_info, rank_one_update
    public :: test_convergence, calculate_rho
    public :: solve_LLS, shift_matrix
    public :: dogleg, more_sorensen, generate_scaling, solve_newton_tensor, aint_tr
    public :: switch_to_quasi_newton
    public :: ERROR
    
contains

  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ S O L V E ****************************!!
  !!******************************************************!!
  !!******************************************************!!

  RECURSIVE SUBROUTINE NLLS_SOLVE( n, m, X,                   & 
                         eval_F, eval_J, eval_HF,   & 
                         params,                    &
                         options, inform, weights )
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( NLLS_inform ), INTENT( OUT ) :: inform
    TYPE( NLLS_options ), INTENT( IN ) :: options
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
    real( wp ), dimension( m ), intent(in), optional :: weights
      
    integer  :: i
    
    type ( NLLS_workspace ) :: w
    
!!$    write(*,*) 'Controls in:'
!!$    write(*,*) 'error = ',options%error
!!$    write(*,*) 'out = ', options%out
!!$    write(*,*) 'print_level = ', options%print_level
!!$    write(*,*) 'maxit = ', options%maxit
!!$    write(*,*) 'model = ', options%model
!!$    write(*,*) 'nlls_method = ', options%nlls_method
!!$    write(*,*) 'lls_solver = ', options%lls_solver
!!$    write(*,*) 'stop_g_absolute = ', options%stop_g_absolute
!!$    write(*,*) 'stop_g_relative = ', options%stop_g_relative     
!!$    write(*,*) 'initial_radius = ', options%initial_radius
!!$    write(*,*) 'maximum_radius = ', options%maximum_radius
!!$    write(*,*) 'eta_successful = ', options%eta_successful
!!$    write(*,*) 'eta_very_successful = ',options%eta_very_successful
!!$    write(*,*) 'eta_too_successful = ',options%eta_too_successful
!!$    write(*,*) 'radius_increase = ',options%radius_increase
!!$    write(*,*) 'radius_reduce = ',options%radius_reduce
!!$    write(*,*) 'radius_reduce_max = ',options%radius_reduce_max
!!$    write(*,*) 'hybrid_switch = ',options%hybrid_switch
!!$    write(*,*) 'subproblem_eig_fact = ',options%subproblem_eig_fact
!!$    write(*,*) 'more_sorensen_maxits = ',options%more_sorensen_maxits
!!$    write(*,*) 'more_sorensen_shift = ',options%more_sorensen_shift
!!$    write(*,*) 'more_sorensen_tiny = ',options%more_sorensen_tiny
!!$    write(*,*) 'more_sorensen_tol = ',options%more_sorensen_tol
!!$    write(*,*) 'hybrid_tol = ', options%hybrid_tol
!!$    write(*,*) 'hybrid_switch_its = ', options%hybrid_switch_its
!!$    write(*,*) 'output_progress_vectors = ',options%output_progress_vectors

    main_loop: do i = 1,options%maxit
       
       if ( present(weights) ) then
          call nlls_iterate(n, m, X,      & 
               w,                         &
               eval_F, eval_J, eval_HF,   & 
               params,                    &
               inform, options, weights)
       else
          call nlls_iterate(n, m, X,      & 
               w,                         &
               eval_F, eval_J, eval_HF,   & 
               params,                    &
               inform, options)
       end if
       ! test the returns to see if we've converged

       if (inform%status < 0) then 
          call nlls_strerror(inform)
          if ( options%print_level > 0 ) then
             write(options%error,'(a,a)') 'ERROR: ', trim(inform%error_message)
          end if
          goto 1000 ! error -- exit
       elseif ((inform%convergence_normf == 1).or.(inform%convergence_normg == 1)) then
          goto 1000 ! converged -- exit
       end if
       
     end do main_loop
    
     ! If we reach here, then we're over maxits     
     if (options%print_level > 0 ) write(options%error,1040) 
     inform%status = ERROR%MAXITS
     goto 1000
    
1000 continue
     call nlls_finalize(w,options)
     return
! Non-executable statements

! print level > 0

1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

   END SUBROUTINE NLLS_SOLVE
  

  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ I T E R A T E ************************!!
  !!******************************************************!!
  !!******************************************************!!

   recursive subroutine nlls_iterate(n, m, X,                   & 
                          w,                         & 
                          eval_F, eval_J, eval_HF,   & 
                          params,                    &
                          inform, options, weights)

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( nlls_inform ), INTENT( INOUT ) :: inform
    TYPE( nlls_options ), INTENT( IN ) :: options
    type( NLLS_workspace ), INTENT( INOUT ) :: w
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
    REAL( wp ), DIMENSION( m ), INTENT( IN ), optional :: weights
      
    integer :: svdstatus = 0
    integer :: i, no_reductions, max_tr_decrease = 100
    real(wp) :: rho, rho_gn, normFnew, md, md_gn, Jmax, JtJdiag
    real(wp) :: FunctionValue, normX, normXnew
    logical :: success 
    logical :: bad_allocate = .false.
    character :: second
    real(wp) :: sum_reg
    integer :: num_successful_steps = 0
    
    ! todo: make max_tr_decrease a control variable

    ! Perform a single iteration of the RAL_NLLS loop
    if (w%first_call == 1) then
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       !! This is the first call...allocate arrays, and get initial !!
       !! function evaluations                                      !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       if ( options%print_level >= 2 ) write( options%out, 2000 ) 
       if ( options%print_level >= 1 ) write( options%out, 1000 )
       ! first, check if n < m
       if (n > m) goto 4070
       if (options%type_of_method == 2) then
          ! now check that if type_of_method = 2, we have an appropriate subproblem solverx
          if (options%nlls_method .ne. 4) goto 4110
       end if
       ! set scalars...
       w%first_call = 0
       w%tr_nu = options%radius_increase
       w%tr_p = 7
       inform%status = 0
       inform%external_return = 0
       ! allocate space for vectors that will be used throughout the algorithm
       
       if (options%setup_workspaces) then 
          call setup_workspaces(w,n,m,options,inform)
          if ( inform%alloc_status > 0) goto 4000
       elseif (.not. w%allocated) then 
          goto 4100
       end if
       
       ! evaluate the residual
       call eval_F(inform%external_return, n, m, X, w%f, params)
       inform%f_eval = inform%f_eval + 1
       if (inform%external_return .ne. 0) goto 4020
       if ( present(weights)) then
          ! set f -> Wf
          w%f(1:m) = weights(1:m)*w%f(1:m)
       end if

       ! and evaluate the jacobian
       call eval_J(inform%external_return, n, m, X, w%J, params)
       inform%g_eval = inform%g_eval + 1
       if (inform%external_return .ne. 0) goto 4010
       if ( present(weights) ) then
          ! set J -> WJ
          do i = 1, n
             w%J( (i-1)*m + 1 : i*m) = weights(1:m)*w%J( (i-1)*m + 1 : i*m)
          end do
       end if
       
       if (options%relative_tr_radius == 1) then 
          ! first, let's get diag(J^TJ)
          Jmax = 0.0
          do i = 1, n
             ! note:: assumes column-storage of J
             JtJdiag = norm2( w%J( (i-1)*m + 1 : i*m ) )
             if (JtJdiag > Jmax) Jmax = JtJdiag
          end do
          w%Delta = options%initial_radius_scale * (Jmax**2)
       else
          w%Delta = options%initial_radius
       end if

       if ( options%calculate_svd_J ) then
          ! calculate the svd of J (if needed)
          call get_svd_J(n,m,w%J,&
               w%smallest_sv(1), w%largest_sv(1), &
               options,inform,svdstatus,w%get_svd_J_ws)
          if ((svdstatus .ne. 0).and.(options%print_level .ge. 3)) then 
             write(options%out,'(a,i0)') 'warning! svdstatus = ', svdstatus
             write( options%out, 3140 ) svdstatus
          end if
       end if

       w%normF = norm2(w%f(1:m))
       if (options%regularization > 0 ) then 
          normX = norm2(X(1:n))
          call update_regularized_normF(w%normF,normX,options)
       end if
       w%normF0 = w%normF

       !    g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g)
       w%g = -w%g
       if (options%regularization > 0) call update_regularized_gradient(w%g,X,normX,options)
       w%normJF = norm2(w%g)
       w%normJF0 = w%normJF
       w%normJFold = w%normJF
       
       ! save some data 
       inform%obj = 0.5 * ( w%normF**2 )
       inform%norm_g = w%normJF
       inform%scaled_g = w%normJF/w%normF

       ! if we need to output vectors of the history of the residual
       ! and gradient, the set the initial values
       if (options%output_progress_vectors) then
          w%resvec(1) = inform%obj
          w%gradvec(1) = inform%norm_g
       end if
       
       ! set the reg_order to that in the options
       w%reg_order = options%reg_order

       !! Select the order of the model to be used..
       select case (options%model)
       case (1) ! first-order
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
          if ( ( options%type_of_method == 2) .and. &
               (options%reg_order .le. zero)) then 
             ! regularization method, use optimal reg
             w%reg_order = two
          end if
       case (2) ! second order
          if ( options%exact_second_derivatives ) then
             if ( present(weights) ) then
                call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
             else
                call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
             end if
             inform%h_eval = inform%h_eval + 1
             if (inform%external_return > 0) goto 4030
             
             if (options%regularization > 0) then 
                call update_regularized_hessian(w%hf,X,n,options)
             end if
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = zero
          end if
          w%use_second_derivatives = .true.
          if ( ( options%type_of_method == 2) .and. & 
               (options%reg_order .le. zero)) then 
             ! regularization method, use optimal reg
             w%reg_order = three
          end if
       case (3) ! hybrid (MNT)
          ! set the tolerance :: make this relative
          w%hybrid_tol = options%hybrid_tol * ( w%normJF/(0.5*(w%normF**2)) )                   
          ! use first-order method initially
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
          if ( (options%type_of_method == 2) .and. & 
               (options%reg_order .le. zero)) then 
             ! regularization method, use optimal reg
             w%reg_order = two
          end if
          if (.not. options%exact_second_derivatives) then 
             ! initialize hf_temp too 
             w%hf_temp(1:n**2) = zero
          end if
          
       case (4) ! tensor model....
          ! get the intitial second derivatives
          w%use_second_derivatives = .true.
          if ( options%exact_second_derivatives ) then
             if ( present(weights) ) then
                call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
             else
                call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
             end if
             inform%h_eval = inform%h_eval + 1
             if (inform%external_return > 0) goto 4030
          else
             goto 4090 ! return an error
          end if
       case default
          goto 4040 ! unsupported model -- return to user
       end select
       
       rho  = -one ! intialize rho as a negative value

       if (options%print_level > 0 ) then 
          write(options%out,1010) w%iter, ' ', w%Delta, rho, inform%obj, &
            inform%norm_g, inform%scaled_g
       end if
       
    end if

    w%iter = w%iter + 1
    inform%iter = w%iter
    

    success = .false.
    no_reductions = 0

    do while (.not. success) ! loop until successful
       no_reductions = no_reductions + 1
       if (no_reductions > max_tr_decrease+1) goto 4050

       if (options%print_level >=1) then
          if (options%model == 4) then
             second = 'T';
          elseif (.not. w%use_second_derivatives) then 
             second = 'N';
          elseif (options%exact_second_derivatives) then 
             second = 'Y';
          else
             second = 'A';
          end if
       end if


       !+++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the step                        !
       !    d                                      !   
       ! that the model thinks we should take next !
       !+++++++++++++++++++++++++++++++++++++++++++!
       w%calculate_step_ws%solve_galahad_ws%reg_order = w%reg_order 
       ! todo: fix this...don't copy over...
       call calculate_step(w%J,w%f,w%hf,w%g,& 
            X,md,md_gn,& 
            n,m,w%use_second_derivatives,w%Delta,eval_HF, params, & 
            num_successful_steps, &
            w%Xnew,w%d,w%normd, & 
            options,inform,& 
            w%calculate_step_ws)
       if (inform%status .ne. 0) goto 4000

       !++++++++++++++++++!
       ! Accept the step? !
       !++++++++++++++++++!
       call eval_F(inform%external_return, n, m, w%Xnew, w%fnew, params)
       inform%f_eval = inform%f_eval + 1

       
       if (inform%external_return .ne. 0) goto 4020
       if ( present(weights) ) then
          ! set f -> Wf
          w%fnew(1:m) = weights(1:m)*w%fnew(1:m)
       end if       
       normFnew = norm2(w%fnew(1:m))
       
       if (options%regularization > 0) then 
          normX = norm2(w%Xnew)
          call update_regularized_normF(normFnew,normX,options)
       end if     
       
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the quantity                                   ! 
       !   rho = 0.5||f||^2 - 0.5||fnew||^2 =   actual_reduction  !
       !         --------------------------   ------------------- !
       !             m_k(0)  - m_k(d)         predicted_reduction !
       !                                                          !
       ! if model is good, rho should be close to one             !
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       call calculate_rho(w%normF,normFnew,md,rho,options)
       if ( (rho > HUGE(wp)) .or. &
            (rho .ne. rho) .or. &
            (rho .le. options%eta_successful) ) then
          ! rho is either very large, NaN, or unsuccessful
          num_successful_steps = 0
          if ( (w%use_second_derivatives) .and.  &
               (options%model == 3) .and. & 
               (no_reductions==1) ) then
             ! recalculate rho based on the approx GN model	
             ! (i.e. the Gauss-Newton model evaluated at the Quasi-Newton step)
             call calculate_rho(w%normF,normFnew,md_gn,rho_gn,options)
             if (rho_gn > options%eta_successful) then
                ! don't trust this model -- switch to GN
                ! (See Dennis, Gay and Walsh (1981), Section 5)
                call switch_to_gauss_newton(w,n,options)
                w%hf_temp(:) = zero
                cycle
             end if
          end if
       else
          ! success!!
          num_successful_steps = num_successful_steps + 1
          success = .true.
       end if
       
       !++++++++++++++++++++++!
       ! Update the TR radius !
       !++++++++++++++++++++++!
       call update_trust_region_radius(rho,options,inform,w)
       if (inform%status .ne. 0) goto 4000
       
       if (.not. success) then
          ! finally, check d makes progress
          if ( options%print_level >= 1 ) then
             write(options%out,1020) w%iter, second, w%Delta, rho
          end if
          if ( norm2(w%d) < epsmch * norm2(w%Xnew) ) goto 4060
       end if
    end do
    ! if we reach here, a successful step has been found
    
    ! update X and f
    X(:) = w%Xnew(:)
    w%f(:) = w%fnew(:)
    
    if (.not. options%exact_second_derivatives) then 
       ! first, let's save some old values...
       ! g_old = -J_k^T r_k
       w%g_old = w%g
       ! g_mixed = -J_k^T r_{k+1}
       call mult_Jt(w%J,n,m,w%fnew,w%g_mixed)
       w%g_mixed = -w%g_mixed
    end if

    ! evaluate J and hf at the new point
    call eval_J(inform%external_return, n, m, X, w%J, params)
    inform%g_eval = inform%g_eval + 1
    if (inform%external_return .ne. 0) goto 4010
    if ( present(weights) ) then
       ! set J -> WJ
       do i = 1, n
          w%J( (i-1)*m + 1 : i*m) = weights(1:m)*w%J( (i-1)*m + 1 : i*m)
       end do
    end if
    
    if ( options%calculate_svd_J ) then
       call get_svd_J(n,m,w%J,&
            w%smallest_sv(w%iter + 1), w%largest_sv(w%iter + 1), &
            options,inform,svdstatus,w%get_svd_J_ws)
       if ((svdstatus .ne. 0).and.(options%print_level >= 3)) then 
          write( options%out, 3140 ) svdstatus
       end if
    end if
    
    ! g = -J^Tf
    call mult_Jt(w%J,n,m,w%f,w%g)
    w%g = -w%g
    if ( options%regularization > 0 ) call update_regularized_gradient(w%g,X,normX,options)
    
    w%normJFold = w%normJF
    w%normF = normFnew
    w%normJF = norm2(w%g)
    
    if (options%model == 3) then
       ! hybrid method -- check if we need second derivatives
       
       if (w%use_second_derivatives) then 
          if (w%normJF > w%normJFold) then 
             ! switch to Gauss-Newton      
             call switch_to_gauss_newton(w,n,options)
          end if
       else
          FunctionValue = 0.5 * (w%normF**2)
          if ( w%normJF/FunctionValue < w%hybrid_tol ) then 
             w%hybrid_count = w%hybrid_count + 1
             if (w%hybrid_count == options%hybrid_switch_its) then
                ! use (Quasi-)Newton
                call switch_to_quasi_newton(w,n,options)
             end if
          else 
             w%hybrid_count = 0
          end if
       end if

       if( .not. w%use_second_derivatives) then
          ! call apply_second_order_info anyway, so that we update the
          ! second order approximation
          if (.not. options%exact_second_derivatives) then
             call rank_one_update(w%hf_temp,w,n)
          end if
       end if

    end if

    if ( w%use_second_derivatives ) then 
       if (present(weights)) then
          call apply_second_order_info(n,m, & 
               X,w, & 
               eval_Hf, &
               params,options,inform, weights)
       else
          call apply_second_order_info(n,m, & 
               X,w, & 
               eval_Hf, &
               params,options,inform)
       end if
       if (inform%external_return .ne. 0) goto 4030
    end if

    ! update the stats 
    inform%obj = 0.5*(w%normF**2)
    inform%norm_g = w%normJF
    inform%scaled_g = w%normJF/w%normF
    if (options%output_progress_vectors) then
       w%resvec (w%iter + 1) = inform%obj
       w%gradvec(w%iter + 1) = inform%norm_g
    end if
    
    if (options%print_level >=1) then
       write(options%out,1010) w%iter, second, w%Delta, rho, inform%obj, &
            inform%norm_g, inform%scaled_g
    end if

    !++++++++++++++++++!
    ! Test convergence !
    !++++++++++++++++++!
    call test_convergence(w%normF,w%normJF,w%normF0,w%normJF0,w%normd,options,inform)
    if (inform%convergence_normf == 1) goto 5000 ! <----converged!!
    if (inform%convergence_normg == 1) goto 5010 ! <----converged!!


! Non-executable statements

! print level > 0

!1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

1000 FORMAT('iter',4x,'2nd order?',2x,'Delta',9x,'rho',11x,'0.5||f||^2',4x,'||J''f||',7x,'||J''f||/||f||')
1010 FORMAT(   i4, 4x,      A,9x,     ES12.4, 2x,ES12.4,2x,ES12.4      ,2x,ES12.4,    2x,ES12.4)
1020 FORMAT(   i4, 4x,      A,9x,     ES12.4, 2x,ES12.4,2x,'unsuccessful step')
! print level > 1
2000 FORMAT(/,'* Running RAL_NLLS *')


! print level > 2
3110 FORMAT('Initial trust region radius taken as ', ES12.4)
3140 FORMAT('Warning: Error when calculating svd, status = ',I0)



! error returns
4000 continue
    ! generic end of algorithm
        ! all (final) exits should pass through here...
    inform%iter = w%iter
!    if (options%output_progress_vectors) then
    if (allocated(w%resvec)) then
       if( allocated(inform%resvec)) deallocate(inform%resvec)
       allocate(inform%resvec(w%iter + 1), stat = inform%alloc_status)
       if (inform%alloc_status > 0) bad_allocate = .true.
       inform%resvec(1:w%iter + 1) = w%resvec(1:w%iter + 1)
    end if
    if (allocated(w%gradvec)) then
       if (allocated(inform%gradvec)) deallocate(inform%gradvec)
       allocate(inform%gradvec(w%iter + 1), stat = inform%alloc_status)
       if (inform%alloc_status > 0) bad_allocate = .true.
       inform%gradvec(1:w%iter + 1) = w%gradvec(1:w%iter + 1)
    end if
    if (options%calculate_svd_J) then
       if (allocated(inform%smallest_sv) ) deallocate(inform%smallest_sv)
       allocate(inform%smallest_sv(w%iter + 1))
       if (inform%alloc_status > 0) bad_allocate = .true.
       if (allocated(inform%largest_sv) ) deallocate(inform%largest_sv)
       allocate(inform%largest_sv(w%iter + 1))
       if (inform%alloc_status > 0) bad_allocate = .true.
    end if

    if (bad_allocate) then 
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = 'nlls_iterate'
    end if

    return

4010 continue
    ! Error in eval_J
    inform%external_name = 'eval_J'
    inform%status = ERROR%EVALUATION
    goto 4000

4020 continue
    ! Error in eval_F
    inform%external_name = 'eval_F'
    inform%status = ERROR%EVALUATION
    goto 4000

4030 continue
    ! Error in eval_HF
    inform%external_name = 'eval_HF'
    inform%status = ERROR%EVALUATION
    goto 4000

4040 continue 
    inform%status = ERROR%UNSUPPORTED_MODEL
    goto 4000

4050 continue 
    ! max tr reductions exceeded
    inform%status = ERROR%MAX_TR_REDUCTIONS
    goto 4000

4060 continue 
    ! x makes no progress
    inform%status = ERROR%X_NO_PROGRESS
    goto 4000

4070 continue
    ! n > m on entry
    inform%status = ERROR%N_GT_M
    goto 4000

4080 continue
    ! bad allocation
    inform%status = ERROR%ALLOCATION
    inform%bad_alloc = 'nlls_iterate'
    goto 4000

4090 continue
    ! no second derivatives in tensor model
    inform%status = ERROR%NO_SECOND_DERIVATIVES
    goto 4000

4100 continue 
    ! workspace error
    inform%status = ERROR%WORKSPACE_ERROR
    goto 4000
    
4110 continue
    ! bad subproblem solver
    inform%status = ERROR%NT_BAD_SUBPROBLEM
    goto 4000

! convergence 
5000 continue
    ! convegence test satisfied
    if (options%print_level >= 2) then
       write(options%out,'(a,i0)') 'RAL_NLLS converged (on ||f|| test) at iteration ', &
            w%iter
    end if
    goto 4000

5010 continue
    if (options%print_level >= 2) then
       write(options%out,'(a,i0)') 'RAL_NLLS converged (on gradient test) at iteration ', &
            w%iter
    end if
    goto 4000

  end subroutine nlls_iterate
  

  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ F I N A L I Z E **********************!!
  !!******************************************************!!
  !!******************************************************!!

  subroutine nlls_finalize(w,options)
    
    type( nlls_workspace ) :: w
    type( nlls_options ) :: options
    
    ! reset all the scalars
    w%first_call = 1
    w%iter = 0
    w%reg_order = 0.0
    w%use_second_derivatives = .false.
    w%hybrid_count = 0 
    w%hybrid_tol = 1.0
    w%tr_nu = 2.0
    w%tr_p = 3
    
    if (options%remove_workspaces) call remove_workspaces(w,options)   

  end subroutine nlls_finalize

  subroutine nlls_strerror(inform)!,error_string)
    type( nlls_inform ), intent(inout) :: inform
    
    if ( inform%status == ERROR%MAXITS ) then
       inform%error_message = 'Maximum number of iterations reached'
    elseif ( inform%status == ERROR%EVALUATION ) then
       write(inform%error_message,'(a,a,a,i0)') & 
            'Error code from user-supplied subroutine ',trim(inform%external_name), & 
            ' passed error = ', inform%external_return
    elseif ( inform%status == ERROR%UNSUPPORTED_MODEL ) then
       inform%error_message = 'Unsupported model passed in options'
    elseif ( inform%status == ERROR%FROM_EXTERNAL ) then
       write(inform%error_message,'(a,a,a,i0)') & 
            'The external subroutine ',trim(inform%external_name), & 
            ' passed error = ', inform%external_return
    elseif ( inform%status == ERROR%UNSUPPORTED_METHOD ) then
       inform%error_message = 'Unsupported nlls_method passed in options'
    elseif ( inform%status == ERROR%ALLOCATION ) then
       write(inform%error_message,'(a,a)') &
            'Bad allocation of memory in ', trim(inform%bad_alloc)
    elseif ( inform%status == ERROR%MAX_TR_REDUCTIONS ) then
       inform%error_message = 'The trust region was reduced the maximum number of times'
    elseif ( inform%status == ERROR%X_NO_PROGRESS ) then
       inform%error_message = 'No progress made in X'
    elseif ( inform%status == ERROR%N_GT_M ) then
       inform%error_message = 'The problem is overdetermined'
    elseif ( inform%status == ERROR%BAD_TR_STRATEGY ) then
       inform%error_message = 'Unsupported tr_update_stategy passed in options'
    elseif ( inform%status == ERROR%FIND_BETA ) then
       inform%error_message = 'Unable to find suitable scalar in findbeta subroutine'
    elseif ( inform%status == ERROR%BAD_SCALING ) then
       inform%error_message = 'Unsupported value of scale passed in options'
    elseif ( inform%status == ERROR%WORKSPACE_ERROR ) then
       inform%error_message = 'Error accessing pre-allocated workspace'
    elseif ( inform%status == ERROR%DOGLEG_MODEL ) then
       inform%error_message = 'Model not supported in dogleg (nlls_method=1)'
    elseif ( inform%status == ERROR%AINT_EIG_IMAG ) then
       inform%error_message = 'All eigenvalues are imaginary (nlls_method=2)'
    elseif ( inform%status == ERROR%AINT_EIG_ODD ) then
       inform%error_message = 'Odd matrix sent to max_eig subroutine (nlls_method=2)'
    elseif ( inform%status == ERROR%MS_MAXITS ) then
       inform%error_message = 'Maximum iterations reached in more_sorensen (nlls_method=3)'
    elseif ( inform%status == ERROR%MS_TOO_MANY_SHIFTS ) then
       inform%error_message = 'Too many shifts taken in more_sorensen (nlls_method=3)'
    elseif ( inform%status == ERROR%MS_NO_PROGRESS ) then
       inform%error_message = 'No progress being made in more_sorensen (nlls_method=3)'
    elseif ( inform%status == ERROR%NO_SECOND_DERIVATIVES ) then
       inform%error_message = 'Exact second derivatives needed for tensor model'
    elseif ( inform%status == ERROR%NT_BAD_SUBPROBLEM ) then
       inform%error_message = 'nlls_method = 4 needed if type_of_method=2'
    else 
       inform%error_message = 'Unknown error number'           
    end if
    
  end subroutine nlls_strerror


! below are the truly internal subroutines...

  RECURSIVE SUBROUTINE calculate_step(J,f,hf,g,X,md,md_gn,n,m,use_second_derivatives, & 
                            Delta,eval_HF,params,num_successful_steps,& 
                            Xnew,d,normd,options,inform,w)

! -------------------------------------------------------
! calculate_step, find the next step in the optimization
! -------------------------------------------------------

    REAL(wp), intent(in) :: J(:), f(:), hf(:)
    REAL(wp), intent(inout) :: g(:)
    REAL(wp), intent(inout) :: Delta
    REAL(wp), intent(out) :: X(:)
    procedure( eval_hf_type ) :: eval_HF       
    class( params_base_type ) :: params  
    integer, intent(in) :: num_successful_steps
    integer, intent(in)  :: n, m
    logical, intent(in) :: use_second_derivatives
    real(wp), intent(out) :: d(:), Xnew(:)
    real(wp), intent(out) :: md, md_gn
    real(wp), intent(out) :: normd
    TYPE( nlls_options ), INTENT( IN ) :: options
    TYPE( nlls_inform ), INTENT( INOUT ) :: inform
    TYPE( calculate_step_work ) :: w
        
    real(wp) :: md_bad
    integer :: i
    logical :: scaling_used = .false.
    real(wp) :: coeff ! coefficient in front of diag(s) if reg. problem being solved
    real(wp) :: sum_reg
    real(wp) :: normx

    if (.not. w%allocated) goto 1010
    
    d(1:n) = zero
    w%scale = one
    
    if ( options%model == 4 ) then
       ! tensor model -- call ral_nlls again
       
       call solve_newton_tensor(J, f, eval_HF, X, n, m, Delta, num_successful_steps,& 
            d, md, params, options, inform, & 
            w%solve_newton_tensor_ws)
       normd = norm2(d(1:n)) ! ||d||_D

       Xnew = X + d
       call evaluate_model(f,J,hf,X,Xnew,d,md_bad,md_gn,m,n,options,inform,w%evaluate_model_ws)
    else 
       ! compute the hessian used in the model 

!       w%scale = zero
       w%extra_scale = zero

       ! Set A = J^T J
       call matmult_inner(J,n,m,w%A)
       ! add any second order information...
       ! so A = J^T J + HF
       w%A = w%A + reshape(hf(1:n**2),[n,n])

       ! and, now, let's add on a reg parameter, if needed
       select case (options%regularization) 
       case (1)
          do i = 1, n             
             w%A(i,i) = w%A(i,i) + options%regularization_term
          end do
          w%extra_scale = options%regularization_term
       case (2)
          if ( .not. allocated(w%xxt) ) allocate (w%xxt(n,n))
          call outer_product(X,n,w%xxt)
          normx = norm2(X(1:n))
          if (normx > epsmch) then
             ! add term from J^TJ
             w%A(1:n,1:n) = w%A(1:n,1:n) + &
                  ( options%regularization_term * options%regularization_power / 2.0 ) * & 
                  normx**(options%regularization_power - 4.0) * w%xxt 
             ! since there's extra terms in the 'real' J, add these to the scaling
             do i = 1, n
                ! add the square of the entries of last row of the 'real' Jacobian
                w%extra_scale(i) = & 
                     (options%regularization_term * options%regularization_power / 2.0 ) * &
                     (normx**(options%regularization_power-4)) * X(i)**2.0 
             end do
          end if
       end select

       ! and let's copy over g to a temp vector v
       ! (it's copied so that we can scale v without losing g)

       w%v(1:n) = -g(1:n)

       ! if scaling needed, do it
       if ( (options%nlls_method == 3) .or. (options%nlls_method == 4) ) then 
          if ( (options%scale .ne. 0) ) then
!             call apply_scaling(J,n,m,w%extra_scale,w%A,w%v, & 
!                  w%apply_scaling_ws,options,inform)
             call generate_scaling(J,w%A,n,m,w%scale,w%extra_scale,& 
                  w%generate_scaling_ws,options,inform)
             scaling_used = .true.
          end if
       end if

       IF (scaling_used) then
          do i = 1,n
             w%v(i) = w%v(i) / w%scale(i)
             w%A(i,:) = w%A(i,:) / w%scale(i)
             w%A(:,i) = w%A(:,i) / w%scale(i)
          end do
       end IF
   

       ! (Gauss-)/(Quasi-)Newton method -- solve as appropriate...

       select case (options%nlls_method)
       case (1) ! Powell's dogleg
          if (options%print_level >= 2) write(options%out,3000) 'dogleg'
          call dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w%dogleg_ws)
        case (2) ! The AINT method
           if (options%print_level >= 2) write(options%out,3000) 'AINT_TR'
           call AINT_TR(J,w%A,f,X,w%v,hf,n,m,Delta,d,normd,options,inform,w%AINT_tr_ws)
        case (3) ! More-Sorensen
           if (options%print_level >= 2) write(options%out,3000) 'More-Sorensen'
           call more_sorensen(w%A,w%v,n,m,Delta,d,normd,options,inform,w%more_sorensen_ws)
        case (4) ! Galahad
           if (options%type_of_method == 1) then
              if (options%print_level >= 2) write(options%out,3000) 'DTRS'
           elseif (options%type_of_method == 2) then 
              if (options%print_level >= 2) write(options%out,3020) 'DRQS'
           end if           
           call solve_galahad(w%A,w%v,n,m,Delta,num_successful_steps, & 
                d,normd,options,inform,w%solve_galahad_ws)
        case default
           inform%status = ERROR%UNSUPPORTED_METHOD
           goto 1000
        end select

        ! reverse the scaling on the step
        if ( (scaling_used) ) then 
           do i = 1, n
              d(i) = d(i) / w%scale(i)
           end do
        end if

        !++++++++++++++++++++++++++++!
        ! Get the value of the model !
        !      md :=   m_k(d)        !
        ! evaluated at the new step  !
        ! and the value of the       !
        ! Gauss-Newton model too     ! 
        !     md := m_k^gn(d)        !
        !++++++++++++++++++++++++++++!        
        Xnew = X + d
        call evaluate_model(f,J,hf,X,Xnew,d,md,md_gn,m,n,options,inform,w%evaluate_model_ws)

     end if



     if (options%print_level >= 2) write(options%out,3010)
         
1000 return
     
return
     
1010 continue 
     inform%status = ERROR%WORKSPACE_ERROR
     return

3000 FORMAT('*** Solving the trust region subproblem using ',A,' ***')
3010 FORMAT('*** Subproblem solution found ***')
3020 FORMAT('*** Solving the regularized subproblem using ',A,' ***')
     

   END SUBROUTINE calculate_step

   subroutine generate_scaling(J,A,n,m,scale,extra_scale,w,options,inform)
     !-------------------------------
     ! generate_scaling
     ! input :: Jacobian matrix, J
     ! ouput :: scaled Hessisan, H, and J^Tf, v.
     !
     ! Calculates a diagonal scaling W, stored in w%diag
     ! updates v(i) -> (1/W_i) * v(i)
     !         A(i,j) -> (1 / (W_i * W_j)) * A(i,j)
     !-------------------------------
     real(wp), intent(in) :: J(*), A(:,:)
     integer, intent(in) :: n,m
     real(wp), intent(inout) :: scale(:), extra_scale(:)
     type( generate_scaling_work ), intent(inout) :: w
     type( nlls_options ), intent(in) :: options
     type( nlls_inform ), intent(inout) :: inform

     integer :: ii, jj
     real(wp) :: Jij, temp

     if (.not. w%allocated) goto 1010

     select case (options%scale)
     case (1,2)
        do ii = 1,n
           temp = extra_scale(ii)
           what_scale: if (options%scale == 1) then
              ! use the scaling present in gsl:
              ! scale by W, W_ii = ||J(i,:)||_2^2
              do jj = 1,m
                 call get_element_of_matrix(J,m,jj,ii,Jij)
                 temp = temp + Jij**2
              end do
           elseif ( options%scale == 2) then 
              ! scale using the (approximate) hessian
              do jj = 1,n
                 temp = temp + A(ii,jj)**2
              end do
           end if what_scale
           trim_scale: if (temp < options%scale_min) then 
              if (options%scale_trim_min) then 
                 temp = options%scale_min
              else
                 temp = one
              end if
           elseif (temp > options%scale_max) then
              if (options%scale_trim_max) then 
                 temp = options%scale_max
              else
                 temp = one 
              end if
           end if trim_scale
           temp = sqrt(temp)
           if (options%scale_require_increase) then
              scale(ii) = max(temp,scale(ii))
           else
              scale(ii) = temp
           end if
        end do
!!$     case (3)
!!$        ! assuming problems are small...do an eigen-decomposition of H
!!$        write(*,*) '***** Warning ********'
!!$        write(*,*) '*    not robust      *'
!!$        write(*,*) '**********************'
!!$        call all_eig_symm(A,n,w%tempvec,w%ev,w%all_eig_symm_ws,inform)
!!$        if (inform%status .ne. 0) goto 1000
!!$        do ii = 1,n
!!$           ! plain version...
!!$           w%diag(n + 1 - ii) = w%tempvec(ii)
!!$        end do
!!$        ! todo : require_increase, test for trimming
     case default
        inform%status = ERROR%BAD_SCALING
        return
     end select
          
     return
     
1000 continue
     ! error in external package
     return

1010 continue
     inform%status = ERROR%WORKSPACE_ERROR
     return
     
   end subroutine generate_scaling

   subroutine switch_to_gauss_newton(w, n, options)
     type (nlls_workspace), intent(inout) :: w
     integer, intent(in) :: n
     type (nlls_options), intent(in) :: options
     
     if (options%print_level .ge. 3) write(options%out,3120) 
     w%use_second_derivatives = .false.
     if ((options%type_of_method == 2) .and. & 
          (options%reg_order .le. zero)) then 
        ! switch to optimal regularization
        w%reg_order = two
     end if
     ! save hf as hf_temp
     w%hf_temp(1:n**2) = w%hf(1:n**2)
     w%hf(1:n**2) = zero
3120 FORMAT('** Switching to Gauss-Newton **')
   end subroutine switch_to_gauss_newton

   subroutine switch_to_quasi_newton(w, n, options)
     type (nlls_workspace), intent(inout) :: w
     integer, intent(in) :: n 
     type (nlls_options), intent(in) :: options
     
     if (options%print_level .ge. 3) write(options%out,3130) 
     w%use_second_derivatives = .true.
     if ((options%type_of_method == 2).and. & 
        (options%reg_order .le. zero)) then
        ! switch to optimal regularization
        w%reg_order = three
     end if
     w%hybrid_count = 0
     ! copy hf from hf_temp
     if (.not. options%exact_second_derivatives) then
        w%hf(1:n**2) = w%hf_temp(1:n**2)
     end if
3130 FORMAT('** Switching to (Quasi-)Newton **')
   end subroutine switch_to_quasi_newton

   SUBROUTINE dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w)
! -----------------------------------------
! dogleg, implement Powell's dogleg method
! -----------------------------------------

     REAL(wp), intent(in) :: J(:), hf(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     TYPE( dogleg_work ) :: w
     
     real(wp) :: alpha, beta

     if (.not. w%allocated ) goto 1010

     !     Jg = J * g
     call mult_J(J,n,m,g,w%Jg)

     alpha = norm2(g)**2 / norm2( w%Jg )**2
       
     w%d_sd = alpha * g;

     ! Solve the linear problem...
     select case (options%model)
     case (1)
        ! linear model...
        call solve_LLS(J,f,n,m,w%d_gn,inform,w%solve_LLS_ws)
        if ( inform%status .ne. 0 ) goto 1000
     case default
        inform%status = ERROR%DOGLEG_MODEL
        return
     end select
     
     if (norm2(w%d_gn) <= Delta) then
        d = w%d_gn
        if (options%print_level >=2) write(options%out,2000)
     else if (norm2( alpha * w%d_sd ) >= Delta) then
        d = (Delta / norm2(w%d_sd) ) * w%d_sd
        if (options%print_level >=2) write(options%out,2010)
     else
        w%d_sd = alpha * w%d_sd
        w%ghat = w%d_gn - w%d_sd
        call findbeta(w%d_sd,w%ghat,Delta,beta,inform)
        if ( inform%status .ne. 0 ) goto 1000
        d = w%d_sd + beta * w%ghat
        if (options%print_level >=2) write(options%out,2020)
     end if

     normd = norm2(d)
     
     return
     
1000 continue 
     ! bad error return from solve_LLS
     return

1010 continue
     inform%status = ERROR%WORKSPACE_ERROR
     return

! Printing commands
2000 FORMAT('Gauss Newton step taken')
2010 FORMAT('Steepest descent step taken')
2020 FORMAT('Dogleg step taken')

   END SUBROUTINE dogleg
     
   SUBROUTINE AINT_tr(J,A,f,X,v,hf,n,m,Delta,d,normd,options,inform,w)
     ! -----------------------------------------
     ! AINT_tr
     ! Solve the trust-region subproblem using 
     ! the method of ADACHI, IWATA, NAKATSUKASA and TAKEDA
     ! -----------------------------------------

     REAL(wp), intent(in) :: J(:), A(:,:), hf(:), f(:), v(:), X(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( AINT_tr_work ) :: w
        
     integer :: keep_p0, i, size_hard(2)
     real(wp) :: obj_p0, obj_p1, obj_p0_gn, obj_p1_gn
     REAL(wp) :: norm_p0, tau, lam, eta

     if ( .not. w%allocated ) goto 1010
     ! todo..
     ! seems wasteful to have a copy of A and B in M0 and M1
     ! use a pointer?

     keep_p0 = 0
     tau = 1e-4
     obj_p0 = HUGE(wp)

     ! The code finds 
     !  min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p||_B \leq Delta
     !
     ! set A and v for the model being considered here...

     ! Set B to I by hand  
     ! todo: make this an option
     w%B = 0
     do i = 1,n
        w%B(i,i) = 1.0
     end do
     
     select case (options%model)
     case (1)
        call solve_spd(A,-v,w%LtL,w%p0,n,inform)
        if (inform%status .ne. 0) goto 1000
     case default
        call solve_general(A,-v,w%p0,n,inform,w%solve_general_ws)
        if (inform%status .ne. 0) goto 1000
     end select
          
     call matrix_norm(w%p0,w%B,norm_p0)
     
     if (norm_p0 < Delta) then
        keep_p0 = 1;
        ! get obj_p0 : the value of the model at p0
        if (options%print_level >=3) write(options%out,2000) 'p0'     
        call evaluate_model(f,J,hf,X,X,w%p0,obj_p0,obj_p0_gn,m,n, & 
             options, inform, w%evaluate_model_ws)
     end if

     w%M0(1:n,1:n) = -w%B
     w%M0(n+1:2*n,1:n) = A
     w%M0(1:n,n+1:2*n) = A
     call outer_product(v,n,w%gtg) ! gtg = Jtg * Jtg^T
     w%M0(n+1:2*n,n+1:2*n) = (-1.0 / Delta**2) * w%gtg

     w%M1 = 0.0
     w%M1(n+1:2*n,1:n) = -w%B
     w%M1(1:n,n+1:2*n) = -w%B
     
     call max_eig(w%M0,w%M1,2*n,lam, w%y, w%y_hardcase, options, inform, w%max_eig_ws)
     if ( inform%status > 0 ) goto 1000

     if (norm2(w%y(1:n)) < tau) then
        ! Hard case
        if ( options%print_level >= 3) write(options%out, 2010)
        ! overwrite H onto M0, and the outer prod onto M1...
        size_hard = shape(w%y_hardcase)
        call matmult_outer( matmul(w%B,w%y_hardcase), size_hard(2), n, w%M1_small)
        w%M0_small = A(:,:) + lam*w%B(:,:) + w%M1_small
        ! solve Hq + g = 0 for q
        select case (options%model) 
        case (1)
           call solve_spd(w%M0_small,-v,w%LtL,w%q,n,inform)
           if (inform%status .ne. 0) goto 1000
        case default
          call solve_general(w%M0_small,-v,w%q,n,inform,w%solve_general_ws)
          if (inform%status .ne. 0) goto 1000
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! (I think..) and inside...fix

        
        ! find max eta st ||q + eta v(:,1)||_B = Delta
        call findbeta(w%q,w%y_hardcase(:,1),Delta,eta,inform)
        if ( inform%status .ne. 0 ) goto 1000

        !!!!!      ^^TODO^^    !!!!!
        ! currently assumes B = I !!
        !!!!       fixme!!      !!!!
        
        w%p1(:) = w%q(:) + eta * w%y_hardcase(:,1)
        
     else 
        select case (options%model)
        case (1)
           call solve_spd(A + lam*w%B,-v,w%LtL,w%p1,n,inform)
           if (inform%status .ne. 0) goto 1000
        case default
           call solve_general(A + lam*w%B,-v,w%p1,n,inform,w%solve_general_ws)
           if (inform%status .ne. 0) goto 1000
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! and inside...fix
     end if
     
     ! get obj_p1 : the value of the model at p1
     if (options%print_level >=3) write(options%out,2000) 'p1'     
     call evaluate_model(f,J,hf,X,X,w%p1,obj_p1,obj_p1_gn,m,n, & 
          options,inform,w%evaluate_model_ws)

     ! what gives the smallest objective: p0 or p1?
     if (obj_p0 < obj_p1) then
        d = w%p0
        if (options%print_level >=2) write(options%out,2030) 'p0'
     else 
        d = w%p1
        if (options%print_level >=2) write(options%out,2030) 'p1'
     end if

     normd = norm2(d)

     return
         
1000 continue 
     ! bad error return from external package
     return

1010 continue
     inform%status = ERROR%WORKSPACE_ERROR
     return    

! print statements   
2000 FORMAT('Evaluating the model at ',A2,':')
2010 FORMAT('Hard case identified')
2030 FORMAT(A2,' chosen as d')
 

   END SUBROUTINE AINT_tr

   subroutine more_sorensen(A,v,n,m,Delta,d,nd,options,inform,w)
     ! -----------------------------------------
     ! more_sorensen
     ! Solve the trust-region subproblem using 
     ! the method of More and Sorensen
     !
     ! Using the implementation as in Algorithm 7.3.6
     ! of Trust Region Methods
     ! 
     ! main output :: d, the soln to the TR subproblem
     ! -----------------------------------------

     REAL(wp), intent(in) :: A(:,:), v(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: nd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ) :: w

     real(wp) :: nq, epsilon
     real(wp) :: sigma, alpha, local_ms_shift, sigma_shift
     integer :: i, no_restarts

     
     ! The code finds 
     !  d = arg min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p|| \leq Delta
     !
     ! set A and v for the model being considered here...

     if (.not. w%allocated) goto 1010
         
     local_ms_shift = options%more_sorensen_shift

     ! d = -A\v
     call solve_spd(A,-v,w%LtL,d,n,inform)
     if (inform%status .eq. 0) then
        ! A is symmetric positive definite....
        sigma = zero
        if (options%print_level >=3) write(options%out,6000)
     else
        ! reset the error calls -- handled in the code....
        inform%status = 0
        inform%external_return = 0
        inform%external_name = REPEAT( ' ', 80 )
        call min_eig_symm(A,n,sigma,w%y1,options,inform,w%min_eig_symm_ws) 
        if (inform%status .ne. 0) goto 1000
        sigma = -(sigma - local_ms_shift)
        if (options%print_level >= 3) write(options%out,6010) sigma
        ! find a shift that makes (A + sigma I) positive definite
        call get_pd_shift(n,A,v,sigma,d,options,inform,w)
        if (inform%status .ne. 0) goto 4000
        if (options%print_level >=3) write(options%out,6020)
     end if
     
     nd = norm2(d)
     
     if (options%print_level >= 2) write(options%out,5000)
     ! now, we're not in the trust region initally, so iterate....
     sigma_shift = zero
     no_restarts = 0
     ! set 'small' in the context of the algorithm
     epsilon = max( options%more_sorensen_tol * Delta, options%more_sorensen_tiny )
     do i = 1, options%more_sorensen_maxits
        if (options%print_level >= 2) write(options%out,5010) i-1, nd, sigma, sigma_shift
                
        if (nd .le. Delta + epsilon) then
           ! we're within the tr radius
           if (options%print_level >= 3) write(options%out,6030)
           if ( abs(sigma) < options%more_sorensen_tiny ) then
              ! we're good....exit
              if (options%print_level >= 3) write(options%out,6040)
              goto 1020
           else if ( abs( nd - Delta ) < epsilon ) then
              ! also good...exit
              if (options%print_level >= 3) write(options%out,6050)
              goto 1020              
           end if
           call findbeta(d,w%y1,Delta,alpha,inform)
           if (inform%status .ne. 0 ) goto 1000  
           d = d + alpha * w%y1
           if (options%print_level >= 3) write(options%out,6060)
           ! also good....exit
           goto 1020
        end if

        w%q = d ! w%q = R'\d
        CALL DTRSM( 'Left', 'Lower', 'No Transpose', 'Non-unit', n, & 
             1, one, w%LtL, n, w%q, n )
        
        nq = norm2(w%q)
        if (options%print_level >= 3) write(options%out,6080) nq
        
        sigma_shift = ( (nd/nq)**2 ) * ( (nd - Delta) / Delta )
        if (abs(sigma_shift) < options%more_sorensen_tiny * abs(sigma) ) then
           if (no_restarts < 1) then 
              ! find a shift that makes (A + sigma I) positive definite
              call get_pd_shift(n,A,v,sigma,d,options,inform,w)
              if (inform%status .ne. 0) goto 4000
              no_restarts = no_restarts + 1
           else
              ! we're not going to make progress...jump out 
              inform%status = ERROR%MS_NO_PROGRESS
              goto 4000
           end if
        else 
           sigma = sigma + sigma_shift
        end if

        call shift_matrix(A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-v,w%LtL,d,n,inform)
        if (inform%status .ne. 0) goto 1000
        
        nd = norm2(d)

     end do
     if (options%print_level >= 2) write(options%out,5010)

     goto 1040
     
1000 continue 
     ! bad error return from external package
     goto 4000
     
1010 continue
     inform%status = ERROR%WORKSPACE_ERROR
     goto 4000

1020 continue
     ! inital point was successful
     if (options%print_level==2) write(options%out,5040)
     goto 4000

1040 continue
     ! maxits reached, not converged
     if (options%print_level >=2) write(options%out,5020)
     inform%status = ERROR%MS_MAXITS
     goto 4000

3000 continue
     ! too many shifts
     inform%status = ERROR%MS_TOO_MANY_SHIFTS
     goto 4000
     
4000 continue
     ! exit the routine
     return 

! Printing statements
! print_level >= 2 
5000 FORMAT('iter',4x,'nd',12x,'sigma',9x,'sigma_shift')
5010 FORMAT(i4,2x,ES12.4,2x,ES12.4,2x,ES12.4)
5020 FORMAT('More-Sorensen failed to converge within max number of iterations')   
5030 FORMAT('More-Sorensen converged at iteration ',i4)
5040 FORMAT('Leaving More-Sorensen')
  
! print_level >= 3 
6000 FORMAT('A is symmetric positive definite')     
6010 FORMAT('Trying a shift of sigma = ',ES12.4)     
6020 FORMAT('A + sigma I is symmetric positive definite') 
6030 FORMAT('We''re within the trust region radius initially')     
6040 FORMAT('Sigma tiny, so exit')  
6050 FORMAT('||d|| = Delta, so exit')   
6060 FORMAT('Return d + alpha*y_1') 
6070 FORMAT('Converged! ||d|| = Delta at iteration ',i4)      
6080 FORMAT('nq = ',ES12.4)

   end subroutine more_sorensen

   subroutine get_pd_shift(n,A,v,sigma,d,options,inform,w)

     !--------------------------------------------------
     ! get_pd_shift
     !
     ! Given an indefinite matrix A, find a shift sigma
     ! such that (A + sigma I) is positive definite
     !--------------------------------------------------
     
     integer, intent(in) :: n 
     real(wp), intent(in) :: A(:,:), v(:)
     real(wp), intent(inout) :: sigma
     real(wp), intent(inout) :: d(:)
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ), intent(inout) :: w

     integer :: no_shifts
     logical :: successful_shift
     
     no_shifts = 0
     successful_shift = .false.
     do while( .not. successful_shift )
        call shift_matrix(A,sigma,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-v,w%LtL,d,n,inform)
        if ( inform%status .ne. 0 ) then
           ! reset the error calls -- handled in the code....
           inform%status = 0
           inform%external_return = 0
           inform%external_name = REPEAT( ' ', 80 )
           no_shifts = no_shifts + 1
           if ( no_shifts == 10 ) goto 3000 ! too many shifts -- exit
           sigma =  sigma + (10**no_shifts) * options%more_sorensen_shift
           if (options%print_level >=3) write(options%out,6010) sigma
        else
           successful_shift = .true.
        end if
     end do

     return

3000 continue
     ! too many shifts
     inform%status = ERROR%MS_TOO_MANY_SHIFTS
     return     

6010 FORMAT('Trying a shift of sigma = ',ES12.4)
     

   end subroutine get_pd_shift
   
   subroutine solve_galahad(A,v,n,m,Delta,num_successful_steps,d,normd,options,inform,w)

     !---------------------------------------------
     ! solve_galahad
     ! Solve the trust-region subproblem using
     ! the DTRS method from Galahad
     ! 
     ! This method needs H to be diagonal, so we need to 
     ! pre-process
     !
     ! main output :: d, the soln to the TR subproblem
     !--------------------------------------------

     REAL(wp), intent(in) :: A(:,:), v(:)
     REAL(wp), intent(inout) :: Delta
     integer, intent(in)  :: n, m
     integer, intent(in) :: num_successful_steps
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D, where D is the scaling
     type( solve_galahad_work ) :: w
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform

     TYPE ( DTRS_CONTROL_TYPE ) :: dtrs_options
     TYPE ( DTRS_inform_type )  :: dtrs_inform
     TYPE ( DRQS_CONTROL_TYPE ) :: drqs_options
     TYPE ( DRQS_inform_type )  :: drqs_inform

!     real(wp), allocatable :: diag(:)
     integer :: ii
     logical :: proceed
     real(wp) :: reg_param

     ! The code finds 
     !  d = arg min_p   w^T p + 0.5 * p^T D p
     !       s.t. ||p|| \leq Delta
     !
     ! where D is diagonal
     !
     ! our probem in naturally in the form
     ! 
     ! d = arg min_p   v^T p + 0.5 * p^T H p
     !       s.t. ||p|| \leq Delta
     !

     if (.not. w%allocated ) goto 1010
     
     ! We have the unprocessed matrices, we need to get an 
     ! eigendecomposition to make A diagonal
     !
     call all_eig_symm(A,n,w%ew,w%ev,w%all_eig_symm_ws,inform)
     if (inform%status .ne. 0) goto 1000

     ! We can now change variables, setting y = Vp, getting
     ! Vd = arg min_(Vx) v^T p + 0.5 * (Vp)^T D (Vp)
     !       s.t.  ||x|| \leq Delta
     ! <=>
     ! Vd = arg min_(Vx) V^Tv^T (Vp) + 0.5 * (Vp)^T D (Vp)
     !       s.t.  ||x|| \leq Delta
     ! <=>

     ! we need to get the transformed vector v
     call mult_Jt(w%ev,n,n,v,w%v_trans)

     ! we've now got the vectors we need, pass to dtrs_solve
      
     do ii = 1,n
        if (abs(w%v_trans(ii)) < epsmch) then
           w%v_trans(ii) = zero
        end if
        if (abs(w%ew(ii)) < epsmch) then
           w%ew(ii) = zero
        end if
     end do
 
    select case (options%type_of_method)
     case (1)
        call dtrs_initialize( dtrs_options, dtrs_inform ) 
        dtrs_options%error = options%error
        dtrs_options%out = options%out
        dtrs_options%print_level = options%print_level - 1
        call dtrs_solve(n, Delta, zero, w%v_trans, w%ew, w%d_trans, & 
                        dtrs_options, dtrs_inform )
        if ( dtrs_inform%status .ne. 0) then
           inform%external_return = dtrs_inform%status
           inform%external_name = 'galahad_dtrs'
           inform%status = ERROR%FROM_EXTERNAL
           goto 1000
        end if
     case(2)
        call drqs_initialize( drqs_options, drqs_inform ) 
        drqs_options%error = options%error
        drqs_options%out = options%out
        drqs_options%print_level = options%print_level - 1
        
        proceed = .false.
        do while (.not. proceed)          
           reg_param = options%base_regularization + 1.0_wp/Delta
           call drqs_solve & 
                (n,w%reg_order, reg_param, zero, w%v_trans, w%ew, w%d_trans, & 
                drqs_options, drqs_inform)
           if ( drqs_inform%status == -7 ) then
              ! drqs_solve has failed because the matrix
              !     J'J + *1/(2*Delta) * I 
              ! is not spd. Fix this by decreasing Delta until it's big enough
              Delta =  Delta / ten
              ! alternatively, we could use a higher order regularization:
              ! call drqs_solve(n,3.0_wp,1.0_wp/Delta, zero, w%v_trans,&
              !                 w%ew, w%d_trans,drqs_options, drqs_inform)
              continue
           elseif ( drqs_inform%status .ne. 0) then
              inform%external_return = drqs_inform%status
              inform%external_name = 'galahad_drqs'
              inform%status = ERROR%FROM_EXTERNAL
              goto 1000
           else
              proceed = .true.
           end if
        end do
     end select
  ! and return the un-transformed vector
     call mult_J(w%ev,n,n,w%d_trans,d)

     normd = norm2(d) ! ||d||_D
     
     return

1000 continue 
     ! bad error return from external package
     return

1010 continue 
     inform%status = ERROR%WORKSPACE_ERROR
     return
     
2000 FORMAT('Regularization order used = ',ES12.4)

   end subroutine solve_galahad


   SUBROUTINE solve_LLS(J,f,n,m,d_gn,inform,w)
       
!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------

       REAL(wp), DIMENSION(:), INTENT(IN) :: J
       REAL(wp), DIMENSION(:), INTENT(IN) :: f
       INTEGER, INTENT(IN) :: n, m
       REAL(wp), DIMENSION(:), INTENT(OUT) :: d_gn
       type(NLLS_inform), INTENT(INOUT) :: inform

       character(1) :: trans = 'N'
       integer :: nrhs = 1, lwork, lda, ldb
       type( solve_LLS_work ) :: w
       
       if (.not. w%allocated) goto 1000
       
       lda = m
       ldb = max(m,n)
       w%temp(1:m) = f(1:m)
       lwork = size(w%work)
       
       w%Jlls(:) = J(:)
       
       call dgels(trans, m, n, nrhs, w%Jlls, lda, w%temp, ldb, w%work, lwork, &
            inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dgels'
          return
       end if

       d_gn = -w%temp(1:n)

       return
       
1000   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return
              
     END SUBROUTINE solve_LLS
     
     SUBROUTINE findbeta(a, b, Delta, beta, inform)

!  -----------------------------------------------------------------
!  findbeta, a subroutine to find the optimal beta such that 
!   || d || = Delta, where d = a + beta * b
!   
!   uses the approach from equation (3.20b), 
!    "Methods for non-linear least squares problems" (2nd edition, 2004)
!    by Madsen, Nielsen and Tingleff      
!  -----------------------------------------------------------------

     real(wp), dimension(:), intent(in) :: a, b 
     real(wp), intent(in) ::  Delta
     real(wp), intent(out) :: beta
     type( nlls_inform ), intent(inout) :: inform
     
     real(wp) :: c, normb2, norma2, discrim, denom
     
     c = dot_product(a,b)

     norma2 = norm2(a)**2
     normb2 = norm2(b)**2

     discrim = c**2 + (normb2)*(Delta**2 - norma2);
     if ( discrim < zero ) then
        inform%status = ERROR%FIND_BETA
        inform%external_name = 'findbeta'
        return
     end if

     if (c .le. 0) then
        beta = (-c + sqrt(discrim) ) / normb2
     else
        beta = (Delta**2 - norma2) / ( c + sqrt(discrim) )
     end if
        

     END SUBROUTINE findbeta
     
     subroutine evaluate_model(f,J,hf,X,Xnew,d,md,md_gn,m,n,options,inform,w)
! --------------------------------------------------
! Input:
! f = f(x_k), J = J(x_k), 
! hf = \sum_{i=1}^m f_i(x_k) \nabla^2 f_i(x_k) (or an approx)
!
! We have a model 
!      m_k(d) = 0.5 f^T f  + d^T J f + 0.5 d^T (J^T J + HF) d
!
! This subroutine evaluates the model at the point d 
! This value is returned as the scalar
!       md :=m_k(d)
! --------------------------------------------------       

       real(wp), intent(in) :: f(:) ! f(x_k)
       real(wp), intent(in) :: d(:) ! direction in which we move
       real(wp), intent(in) :: J(:) ! J(x_k) (by columns)
       real(wp), intent(in) :: hf(:)! (approx to) \sum_{i=1}^m f_i(x_k) \nabla^2 f_i(x_k)
       real(wp), intent(in) :: X(:) ! original step
       real(wp), intent(in) :: Xnew(:) ! updated step
       integer, intent(in) :: m,n
       real(wp), intent(out) :: md  ! m_k(d)
       real(wp), intent(out) :: md_gn
       TYPE( nlls_options ), INTENT( IN ) :: options
       TYPE( nlls_inform ), INTENT( INOUT ) :: inform
       type( evaluate_model_work ) :: w

       real(wp) :: xtx, xtd, dtd, normx, p, sigma
       
       if (.not. w%allocated ) goto 2000
       md = zero
       md_gn = zero

       !Jd = J*d
       call mult_J(J,n,m,d,w%Jd)
       
       md_gn = 0.5 * norm2(f(1:m) + w%Jd(1:m))**2
       
       ! if we are solving a regularized problem, update terms
       p = options%regularization_power
       sigma = options%regularization_term
       select case (options%regularization) 
       case (1)
          md_gn = md_gn + & 
               0.5 * sigma * norm2(Xnew(1:n))**2
       case (2)
          normx = norm2(X(1:n))
          xtx = normx**2
          xtd = dot_product(X(1:n),d(1:n))
          md_gn = md_gn + & 
               sigma * ( one/p * (normx**p) + & 
               (normx**(p-2)) * xtd + & 
               (p/4.0_wp) * (normx**(p-4)) * (xtd**2) ) 
       end select
       
       
       select case (options%model)
       case (1) ! first-order (no Hessian)
          md = md_gn
          continue
       case (4) ! tensor model 
          ! nothing to do here...
          continue
       case default
          ! these have a dynamic H -- recalculate
          ! H = J^T J + HF, HF is (an approx?) to the Hessian
          call mult_J(hf,n,n,d,w%Hd)
          md = md_gn + 0.5 * dot_product(d(1:n),w%Hd(1:n))
          ! regularized newton terms taken care of already in apply_second_order_info
       end select
       if (options%print_level >= 3) write(options%out,1000) md

       return

1000   FORMAT('Model evauated successfully.  m_k(d) = ',ES12.4)
2000   continue 
       inform%status = ERROR%WORKSPACE_ERROR
       return


     end subroutine evaluate_model

     subroutine calculate_rho(normf,normfnew,md,rho,options)
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the quantity                                   ! 
       !   rho = 0.5||f||^2 - 0.5||fnew||^2 =   actual_reduction  !
       !         --------------------------   ------------------- !
       !             m_k(0)  - m_k(d)         predicted_reduction !
       !                                                          !
       ! if model is good, rho should be close to one             !
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       real(wp), intent(in)  :: normf    ! ||f(x_k)|| at 
       real(wp), intent(in)  :: normfnew ! ||f(x_k + d)||
       real(wp), intent(in)  :: md       !    m_k(d)
       real(wp), intent(out) :: rho      ! act_red / pred_red (close to 1 == good)
       TYPE( nlls_options ), INTENT( IN ) :: options

       real(wp) :: actual_reduction, predicted_reduction
       real(wp) :: tol
       
       actual_reduction = ( 0.5 * (normf**2) ) - ( 0.5 * (normfnew**2) )
       predicted_reduction = ( ( 0.5 * (normf**2) ) - md )

       tol = 10 * epsmch
       
       if ( abs(actual_reduction) < tol ) then 
          rho = one
       else if (abs( predicted_reduction ) < tol ) then 
          rho = one
       else
          rho = actual_reduction / predicted_reduction
       end if

       if (options%print_level >= 3) write(options%out,1000) actual_reduction
       if (options%print_level >= 3) write(options%out,1010) predicted_reduction
       if (options%print_level >= 3) write(options%out,1020) rho
       
1000   FORMAT('Actual reduction (in cost function) = ', ES12.4)
1010   FORMAT('Predicted reduction (in model) = ', ES12.4)
1020   FORMAT('rho returned = ', ES12.4)       

     end subroutine calculate_rho

     subroutine apply_second_order_info(n,m,&
          X,w, & !f,hf,
          eval_Hf,&
          !          d, y, y_sharp, & 
          params,options,inform,weights)
       integer, intent(in)  :: n, m 
       real(wp), intent(in) :: X(:)
       type( NLLS_workspace ), intent(inout) :: w
       procedure( eval_hf_type ) :: eval_Hf
       class( params_base_type ) :: params
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       real(wp), intent(in), optional :: weights(:)
       
       if (options%exact_second_derivatives) then
          if ( present(weights) ) then
             call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
          else
             call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
          end if
          inform%h_eval = inform%h_eval + 1
       else
          ! use the rank-one approximation...
          call rank_one_update(w%hf,w,n)                      
       end if
       
       ! update the hessian here if we're solving a regularized problem
       if (options%regularization > 0) then
          call update_regularized_hessian(w%hf,X,n,options)
       end if
       
     end subroutine apply_second_order_info

     subroutine update_regularized_normF(normF,normX,options)
       real(wp), intent(inout) :: normF
       real(wp), intent(in) :: normX
       type( nlls_options ), intent(in) :: options
       
       select case(options%regularization)
       case (1)
          normF = sqrt(normF**2 + & 
               options%regularization_term *  & 
               normX**2 )
       case (2)
          normF = sqrt(normF**2 + & 
               ( 2 * options%regularization_term / options%regularization_power )  *  & 
               normX**options%regularization_power )
       end select
       
     end subroutine update_regularized_normF

     subroutine update_regularized_gradient(g,X,normX,options)
       real(wp), intent(inout) :: g(:)
       real(wp), intent(in) :: X(:), normX
       type( nlls_options ), intent(in) :: options
             
       select case(options%regularization)
       case (1)
          g = g - options%regularization_term * X
       case (2)
          g = g - options%regularization_term *  & 
               (normX**(options%regularization_power - 2.0_wp)) * X 
       end select

     end subroutine update_regularized_gradient

     subroutine update_regularized_hessian(hf,X,n,options)
       real(wp), intent(inout) :: hf(:)
       real(wp), intent(in) :: X(:)
       integer, intent(in) :: n
       type( nlls_options), intent(in) :: options

       integer :: ii, jj        
       real(wp) :: hf_local, p, sigma, normx

       ! given a Hessian, hf, update as required if solving a regularized problem
       if (options%regularization == 2 ) then 
          p = options%regularization_power
          sigma = options%regularization_term
          normx = norm2(X(1:n))
          do ii = 1,n
             do jj = 1,n
                hf_local = x(ii)*x(jj)
                if (ii == jj) hf_local = hf_local + normx**2
                hf_local = sigma * normx**(p - 4.0) * hf_local                
                hf( (ii-1)*n + jj) = hf( (ii-1)*n + jj) + & 
                     hf_local
             end do
          end do
       end if
      
     end subroutine update_regularized_hessian

     subroutine rank_one_update(hf,w,n)

       real(wp), intent(inout) :: hf(:)
       type( NLLS_workspace ), intent(inout) :: w
       integer, intent(in) :: n
      
       real(wp) :: yts, alpha, dSks

       w%y       = w%g_old   - w%g
       w%y_sharp = w%g_mixed - w%g

       yts = dot_product(w%d,w%y)
       if ( abs(yts) < 10 * epsmch ) then
          ! safeguard: skip this update
          return
       end if

       call mult_J(hf,n,n,w%d,w%Sks) ! hfs = S_k * d

       w%ysharpSks = w%y_sharp - w%Sks

       ! now, let's scale hd (Nocedal and Wright, Section 10.2)
       dSks = abs(dot_product(w%d,w%Sks))
       if ( abs(dSks) < 10.0 * epsmch ) then
          ! check this first to avoid possible overflow
          alpha = one
       else
          alpha = abs(dot_product(w%d,w%y_sharp))/ dSks
          alpha = min(one,alpha)
       end if
       hf(:)  = alpha * hf(:)

       ! update S_k (again, as in N&W, Section 10.2)

       ! hf = hf + (1/yts) (y# - Sk d)^T y:
       alpha = 1/yts
       call dGER(n,n,alpha,w%ysharpSks,1,w%y,1,hf,n)
       ! hf = hf + (1/yts) y^T (y# - Sk d):
       call dGER(n,n,alpha,w%y,1,w%ysharpSks,1,hf,n)
       ! hf = hf - ((y# - Sk d)^T d)/((yts)**2)) * y y^T
       alpha = -dot_product(w%ysharpSks,w%d)/(yts**2)
       call dGER(n,n,alpha,w%y,1,w%y,1,hf,n)

     end subroutine rank_one_update


     subroutine update_trust_region_radius(rho,options,inform,w)

       real(wp), intent(inout) :: rho ! ratio of actual to predicted reduction
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_workspace ), intent(inout) :: w

       select case(options%tr_update_strategy)
       case(1) ! default, step-function
          if (rho < options%eta_success_but_reduce) then
             ! unsuccessful....reduce Delta
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             if (options%print_level > 2) write(options%out,3010) w%Delta     
          else if (rho < options%eta_very_successful) then 
             ! doing ok...retain status quo
             if (options%print_level > 2) write(options%out,3020) w%Delta 
          else if (rho < options%eta_too_successful ) then
             ! more than very successful -- increase delta
             select case(options%type_of_method)
             case(1)
                w%Delta = min(options%maximum_radius, &
                     options%radius_increase * w%normd )
                ! if we have a trust region method, then we 
                ! increase based on ||d||, not on Delta, as there's 
                ! no point increasing the radius if we're within the 
                ! trust region
             case(2)
                w%Delta = min(options%maximum_radius, &
                     options%radius_increase * w%Delta )
                ! increase based on Delta to ensure the 
                ! regularized case works too
             end select
             if (options%print_level > 2) write(options%out,3030) w%Delta
          else if (rho >= options%eta_too_successful) then
             ! too successful....accept step, but don't change w%Delta
             if (options%print_level > 2) write(options%out,3040) w%Delta 
          else
             ! just incase (NaNs and the like...)
             if (options%print_level > 2) write(options%out,3050) w%Delta 
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             rho = -one ! set to be negative, so that the logic works....
          end if
       case(2) ! Continuous method
          ! Based on that proposed by Hans Bruun Nielsen, TR IMM-REP-1999-05
          ! http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf
          if (rho >= options%eta_too_successful) then
             ! too successful....accept step, but don't change w%Delta
             if (options%print_level > 2) write(options%out,3040) w%Delta 
          else if (rho > options%eta_successful) then 
             w%Delta = w%Delta * min(options%radius_increase, &
                  max(options%radius_reduce, & 
                  1 - ( (options%radius_increase - 1) * ((1 - 2*rho)**w%tr_p)) ))
             w%tr_nu = options%radius_reduce
             if (options%print_level > 2) write(options%out,3060) w%Delta 
          else if ( rho <= options%eta_successful ) then 
             w%Delta = w%Delta * w%tr_nu
             w%tr_nu =  w%tr_nu * 0.5_wp
             if (options%print_level > 2) write(options%out,3010) w%Delta
          else
             ! just incase (NaNs and the like...)
             if (options%print_level > 2) write(options%out,3050) w%Delta 
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta 
             rho = -one ! set to be negative, so that the logic works....
          end if
       case default
          inform%status = ERROR%BAD_TR_STRATEGY
          return          
       end select

       return

       ! print statements

3010   FORMAT('Unsuccessful step -- decreasing Delta to', ES12.4)      
3020   FORMAT('Successful step -- Delta staying at', ES12.4)     
3030   FORMAT('Very successful step -- increasing Delta to', ES12.4)
3040   FORMAT('Step too successful -- Delta staying at', ES12.4) 
3050   FORMAT('NaN encountered -- reduced Delta to', ES12.4)   
3060   FORMAT('Changing Delta to ', ES12.4)


     end subroutine update_trust_region_radius

     subroutine test_convergence(normF,normJF,normF0,normJF0,normd,options,inform)

       real(wp), intent(in) :: normF, normJf, normF0, normJF0, normd
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform

       if ( normF <= max(options%stop_f_absolute, &
            options%stop_f_relative * normF0) ) then
          inform%convergence_normf = 1
          return
       end if

       if ( (normJF/normF) <= max(options%stop_g_absolute, &
            options%stop_g_relative * (normJF0/normF0)) ) then
          inform%convergence_normg = 1
       end if

       if ( normd < options%stop_s ) then
          inform%convergence_norms = 1
       end if

       return

     end subroutine test_convergence

     subroutine mult_J(J,n,m,x,Jx)
       real(wp), intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       real(wp), intent(out) :: Jx(*)

       real(wp) :: alpha, beta

       Jx(1:m) = 1.0
       alpha = 1.0
       beta  = 0.0

       call dgemv('N',m,n,alpha,J,m,x,1,beta,Jx,1)

     end subroutine mult_J

     subroutine mult_Jt(J,n,m,x,Jtx)
       double precision, intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       double precision, intent(out) :: Jtx(*)

       double precision :: alpha, beta

       Jtx(1:n) = one
       alpha = one
       beta  = zero

       call dgemv('T',m,n,alpha,J,m,x,1,beta,Jtx,1)

     end subroutine mult_Jt

     subroutine get_element_of_matrix(J,m,ii,jj,Jij)
       real(wp), intent(in) :: J(*)
       integer, intent(in) :: m
       integer, intent(in) :: ii,jj
       real(wp), intent(out) :: Jij

       ! return the (ii,jj)th entry of a matrix 

       ! J held by columns....
       Jij = J(ii + (jj-1)*m)

     end subroutine get_element_of_matrix

     subroutine solve_spd(A,b,LtL,x,n,inform)
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: LtL(:,:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform), intent(inout) :: inform

       ! A wrapper for the lapack subroutine dposv.f
       ! get workspace for the factors....
       LtL(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dposv('L', n, 1, LtL, n, x, n, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dposv'
          return
       end if

     end subroutine solve_spd

     subroutine solve_general(A,b,x,n,inform,w)
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform ), intent(inout) :: inform
       type( solve_general_work ) :: w

       ! A wrapper for the lapack subroutine dposv.f
       ! NOTE: A would be destroyed
       
       if (.not. w%allocated) goto 1000
       
       w%A(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dgesv( n, 1, w%A, n, w%ipiv, x, n, inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dgesv'
          return
       end if

       return

       1000   continue ! workspace error
       inform%status = ERROR%WORKSPACE_ERROR
       return


     end subroutine solve_general

     subroutine matrix_norm(x,A,norm_A_x)
       REAL(wp), intent(in) :: A(:,:), x(:)
       REAL(wp), intent(out) :: norm_A_x

       ! Calculates norm_A_x = ||x||_A = sqrt(x'*A*x)

       norm_A_x = sqrt(dot_product(x,matmul(A,x)))

     end subroutine matrix_norm

     subroutine matmult_inner(J,n,m,A)

       integer, intent(in) :: n,m 
       real(wp), intent(in) :: J(*)
       real(wp), intent(out) :: A(n,n)

       ! Takes an m x n matrix J and forms the 
       ! n x n matrix A given by
       ! A = J' * J

       call dgemm('T','N',n, n, m, one,&
            J, m, J, m, & 
            zero, A, n)


     end subroutine matmult_inner

     subroutine matmult_outer(J,n,m,A)

       integer, intent(in) :: n,m 
       real(wp), intent(in) :: J(*)
       real(wp), intent(out) :: A(m,m)

       ! Takes an m x n matrix J and forms the 
       ! m x m matrix A given by
       ! A = J * J'

       call dgemm('N','T',m, m, n, one,&
            J, m, J, m, & 
            zero, A, m)


     end subroutine matmult_outer

     subroutine outer_product(x,n,xxt)

       real(wp), intent(in) :: x(:)
       integer, intent(in) :: n
       real(wp), intent(out) :: xxt(:,:)

       ! Takes an n vector x and forms the 
       ! n x n matrix xtx given by
       ! xxt = x * x'

       xxt(1:n,1:n) = zero
       call dger(n, n, one, x, 1, x, 1, xxt, n)

     end subroutine outer_product

     subroutine all_eig_symm(A,n,ew,ev,w,inform)
       ! calculate all the eigenvalues of A (symmetric)

       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew(:), ev(:,:)
       type( all_eig_symm_work ) :: w
       type( nlls_inform ), intent(inout) :: inform

       integer :: lwork

       if (.not. w%allocated) goto 1000

       ! copy the matrix A into the eigenvector array
       ev(1:n,1:n) = A(1:n,1:n)

       lwork = size(w%work)
       ! call dsyev --> all eigs of a symmetric matrix

       call dsyev('V', & ! both ew's and ev's 
            'U', & ! upper triangle of A
            n, ev, n, & ! data about A
            ew, w%work, lwork, & 
            inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dsyev'
          return
       end if

       return

1000   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return
       
     end subroutine all_eig_symm

     subroutine min_eig_symm(A,n,ew,ev,options,inform,w)
       ! calculate the leftmost eigenvalue of A

       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew, ev(:)
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_options ), INTENT( IN ) :: options
       type( min_eig_symm_work ) :: w

       real(wp) :: tol, dlamch
       integer :: lwork, eigsout, minindex(1)

       if ( .not. w%allocated ) goto 1000

       tol = 2*dlamch('S')!1e-15

       w%A(1:n,1:n) = A(1:n,1:n) ! copy A, as workspace for dsyev(x)
       ! note that dsyevx (but not dsyev) only destroys the lower (or upper) part of A
       ! so we could possibly reduce memory use here...leaving for 
       ! ease of understanding for now.

       lwork = size(w%work)
       if ( options%subproblem_eig_fact ) then
          ! call dsyev --> all eigs of a symmetric matrix
          call dsyev('V', & ! both ew's and ev's 
               'U', & ! upper triangle of A
               n, w%A, n, & ! data about A
               w%ew, w%work, lwork, & 
               inform%external_return)
          if (inform%external_return .ne. 0) then
             inform%status = ERROR%FROM_EXTERNAL
             inform%external_name = 'lapack_dsyev'
          end if
          minindex = minloc(w%ew)
          ew = w%ew(minindex(1))
          ev = w%A(1:n,minindex(1))
       else
          ! call dsyevx --> selected eigs of a symmetric matrix
          call dsyevx( 'V',& ! get both ew's and ev's
               'I',& ! just the numbered eigenvalues
               'U',& ! upper triangle of A
               n, w%A, n, & 
               1.0, 1.0, & ! not used for RANGE = 'I'
               1, 1, & ! only find the first eigenpair
               tol, & ! abstol for the eigensolver
               eigsout, & ! total number of eigs found
               ew, ev, & ! the eigenvalue and eigenvector
               n, & ! ldz (the eigenvector array)
               w%work, lwork, w%iwork, &  ! workspace
               w%ifail, & ! array containing indicies of non-converging ews
               inform%external_return)
          if (inform%external_return .ne. 0) then
             inform%status = ERROR%FROM_EXTERNAL
             inform%external_name = 'lapack_dsyevx'
          end if
       end if

       return

1000   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return

     end subroutine min_eig_symm

     subroutine max_eig(A,B,n,ew,ev,nullevs,options,inform,w)

       real(wp), intent(inout) :: A(:,:), B(:,:)
       integer, intent(in) :: n 
       real(wp), intent(out) :: ew, ev(:)
       real(wp), intent(inout), allocatable :: nullevs(:,:)
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( max_eig_work ) :: w

       integer :: lwork, maxindex(1), no_null, halfn
       real(wp):: tau
       integer :: i 

       if (.not. w%allocated) goto 2010

       ! Find the max eigenvalue/vector of the generalized eigenproblem
       !     A * y = lam * B * y
       ! further, if ||y(1:n/2)|| \approx 0, find and return the 
       ! eigenvectors y(n/2+1:n) associated with this

       ! check that n is even (important for hard case -- see below)
       if (modulo(n,2).ne.0) goto 1010

       halfn = n/2
       lwork = size(w%work)
       call dggev('N', & ! No left eigenvectors
            'V', &! Yes right eigenvectors
            n, A, n, B, n, &
            w%alphaR, w%alphaI, w%beta, & ! eigenvalue data
            w%vr, n, & ! not referenced
            w%vr, n, & ! right eigenvectors
            w%work, lwork, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = ERROR%FROM_EXTERNAL
          inform%external_name = 'lapack_dggev'
          return
       end if

       ! now find the rightmost real eigenvalue
       w%vecisreal = .true.
       where ( abs(w%alphaI) > 1e-8 ) w%vecisreal = .false.

       w%ew_array(:) = w%alphaR(:)/w%beta(:)
       maxindex = maxloc(w%ew_array,w%vecisreal)
       if (maxindex(1) == 0) goto 1000

       tau = 1e-4 ! todo -- pass this through from above...
       ! note n/2 always even -- validated by test on entry
       if (norm2( w%vr(1:halfn,maxindex(1)) ) < tau) then 
          ! hard case
          ! let's find which ev's are null...
          w%nullindex = 0
          no_null = 0
          do i = 1,n
             if (norm2( w%vr(1:halfn,i)) < 1e-4 ) then
                no_null = no_null + 1 
                w%nullindex(no_null) = i
             end if
          end do
!          allocate(nullevs(halfn,no_null))
          if (no_null > size(nullevs,2)) then
             ! increase the size of the allocated array only if we need to
             if(allocated(nullevs)) deallocate( nullevs )
             allocate( nullevs(halfn,no_null) , stat = inform%alloc_status)
             if (inform%alloc_status > 0) goto 2000
          end if 
          nullevs(1:halfn,1:no_null) = w%vr(halfn+1 : n,w%nullindex(1:no_null))
       end if

       ew = w%alphaR(maxindex(1))/w%beta(maxindex(1))
       ev(:) = w%vr(:,maxindex(1))

       return 

1000   continue 
       inform%status = ERROR%AINT_EIG_IMAG ! Eigs imaginary error
       return

1010   continue
       inform%status = ERROR%AINT_EIG_ODD
       return
       
2000   continue
       inform%status = ERROR%ALLOCATION
       inform%bad_alloc = "max_eig"
       return

2010   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return

     end subroutine max_eig

     subroutine shift_matrix(A,sigma,AplusSigma,n)

       real(wp), intent(in)  :: A(:,:), sigma
       real(wp), intent(out) :: AplusSigma(:,:)
       integer, intent(in) :: n 

       integer :: i 
       ! calculate AplusSigma = A + sigma * I

       AplusSigma(:,:) = A(:,:)
       do i = 1,n
          AplusSigma(i,i) = AplusSigma(i,i) + sigma
       end do

     end subroutine shift_matrix

     subroutine get_svd_J(n,m,J,s1,sn,options,inform,status,w)
       integer, intent(in) :: n,m 
       real(wp), intent(in) :: J(:)
       real(wp), intent(out) :: s1, sn
       type( nlls_options ) :: options
       type( nlls_inform ), intent(inout) :: inform
       integer, intent(out) :: status
       type( get_svd_J_work ) :: w

       !  Given an (m x n)  matrix J held by columns as a vector,
       !  this routine returns the largest and smallest singular values
       !  of J.

       character :: jobu(1), jobvt(1)
       integer :: lwork

       if (.not. w%allocated) goto 1000

       w%Jcopy(:) = J(:)

       jobu  = 'N' ! calculate no left singular vectors
       jobvt = 'N' ! calculate no right singular vectors

       lwork = size(w%work)

       call dgesvd( JOBU, JOBVT, n, m, w%Jcopy, n, w%S, w%S, 1, w%S, 1, & 
            w%work, lwork, status )
       if ( (status .ne. 0) .and. (options%print_level > 3) ) then 
          write(options%error,'(a,i0)') 'Error when calculating svd, dgesvd returned', &
               status
          s1 = -1.0
          sn = -1.0
          ! allow to continue, but warn user and return zero singular values
       else
          s1 = w%S(1)
          sn = w%S(n)
          if (options%print_level > 2) then 
             write(options%out,'(a,es12.4,a,es12.4)') 's1 = ', s1, '    sn = ', sn
             write(options%out,'(a,es12.4)') 'k(J) = ', s1/sn
          end if
       end if

       return
       
1000   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return

     end subroutine get_svd_J


     ! routines needed for the Newton tensor model
     
     subroutine solve_newton_tensor(J, f, eval_HF, X, n, m, Delta, & 
                                    num_successful_steps, & 
                                    d, md, params, options, inform, & 
                                    w)
       
       integer, intent(in)   :: n,m 
       real(wp), intent(in)  :: f(:), J(:)
       real(wp) , intent(in) :: X(:), Delta
       integer, intent(in) :: num_successful_steps
       real(wp), intent(out) :: d(:)
       real(wp), intent(out) :: md
       procedure( eval_hf_type ) :: eval_HF
       class( params_base_type ) :: params              
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( solve_newton_tensor_work ) :: w
              
       type( nlls_inform ) :: tensor_inform
       integer :: i, m_in     
    
       ! We need to solve the problem 
       !   min 1/2 \sum_{i=1}^m t_{ik}^2(s) + 1/p \sigma_k ||s||^p_p
       ! where 
       !   t_{ik}(s) := r_i(x_k) + s' g_i(x_k) + 1/2 s' B_ik s 
       ! and where B_ik is a symmetric approx to Hi(x), the Hessian of r_i(x).
       ! 
       ! To do this, we call ral_nlls recursively.

       ! First, we need to set up the eval_r/J/Hf functions needed here.
       
       if (.not. w%allocated) goto 1000

       ! save the residual to params
       w%tparams%f(1:m) = f(1:m)

       w%tparams%Delta = Delta

       ! save the Jacobian to params
       w%tparams%J(1:n*m) = J(1:n*m)

       ! now, let's get all the Hi's...
       do i = 1,m
          call get_Hi(n, m, X, params, i, w%tparams%Hi(:,:,i), eval_HF, inform)
       end do

       d(1:n) = zero

       ! send to ral_nlls to solve the subproblem recursively
       if (options%print_level > 0) write(options%out,"(80('*'))")
       select case (options%inner_method)
       case (1) ! send in a base regularization parameter
          w%tensor_options%base_regularization = 1.0_wp/Delta
       case (2)
          ! do nothing!
          w%tparams%p = 3
          w%tparams%extra = 2
          w%m_in = m + 1
          d(1:n) = 1e-12 ! Hessian not defined at 0 if p /= 2, so set 'small'          
       case (3)
          w%tensor_options%regularization_term = 1.0_wp / Delta
       case (4)
          d(1:n) = 1e-12 ! Hessian not defined at 0 if p /= 2, so set 'small'
          w%tensor_options%regularization_term = 1.0_wp / Delta
       end select
       do i = 1, w%tensor_options%maxit
          call nlls_iterate(n,w%m_in,d, & 
               inner_workspace, & 
               evaltensor_f, evaltensor_J, evaltensor_HF, &
               w%tparams, & 
               tensor_inform, w%tensor_options )
          if (tensor_inform%status < 0) then
             ! there's an error : exit
             exit
          elseif ( (tensor_inform%convergence_normf == 1) & 
               .or.(tensor_inform%convergence_normg == 1)) then
             ! we've converged!
             exit
          end if
       end do
       call nlls_finalize(inner_workspace,w%tensor_options)
       if (tensor_inform%status .ne. 0) then
          write(options%out,'(A)') '**** inner iteration did not converge ****'
       end if
       inform%inner_iter = inform%inner_iter + tensor_inform%iter
       if (options%print_level > 0) then
          write(options%out,'(a,i0)') 'Total inner iterations = ', inform%inner_iter
          write(options%out,"(80('*'))")
       end if
        
        ! now we need to evaluate the model at the new point
        w%tparams%extra = 0
        call evaltensor_f(inform%external_return, n, m, d, &
             w%model_tensor, w%tparams)
        md = 0.5 * norm2( w%model_tensor(1:m) )**2
        ! + 0.5 * (1.0/Delta) * (norm2(d(1:n))**2)


       return
       
1000   continue
       inform%status = ERROR%WORKSPACE_ERROR
       return
       
     end subroutine solve_newton_tensor

     subroutine evaltensor_f(status, n, m, s, f, params)
       
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(out)   :: f
       class( params_base_type ), intent(in) :: params

       real(wp) :: t_ik
       integer :: ii, jj, kk
       
       ! note:: tenJ is a global (but private to this module) derived type 
       !        with components Hs and Js, which are real arrays 

       ! The function we need to minimize is 
       !  \sum_{i=1}^m t_ik(s) = 1/2 \sum_{i=1}^m (r_i(x_l) + s' g_i(x_k) + 1/2 s' B_ik s)^2

       select type(params)
       type is(tensor_params_type)
          ! note: params%m contains 'm' from the original problem
          ! if we're passing in the reg. factor via the function/Jacobian, then 
          ! 'm' here is m+n from the original problem
          f(1:params%m) = params%f(1:params%m)  ! f_tensor = r_nlls originally
          
          call mult_J(params%J(1:n*params%m),n,params%m,s,tenJ%Js)
          f(1:params%m) = f(1:params%m) + tenJ%Js(1:params%m) ! f_tensor = r + J s
          do ii = 1,params%m
             t_ik = params%f(ii)
             t_ik = t_ik + dot_product(s(1:n),params%J(ii : n*params%m : params%m))
             tenJ%Hs(1:n) = zero
             do jj = 1, n
                ! get Hs = H_i * s
                do kk = 1,n
                   tenJ%Hs(jj) = tenJ%Hs(jj) +  params%Hi(jj,kk,ii) * s(kk)
                end do
             end do
             ! f_tensor_i = r_i + J_is + 0.5 * s'H_is
             f(ii) = f(ii) + 0.5 * dot_product(s(1:n),tenJ%Hs(1:n))
          end do
          if (params%extra == 1) then 
             ! we're passing in the regularization via the function/Jacobian
             do ii = 1, n
                f(params%m + ii) = sqrt(1.0/params%Delta) * s(ii)
             end do
          elseif (params%extra == 2) then 
             f(params%m + 1) = sqrt(2.0/(params%Delta * params%p)) * & 
                               (norm2(s(1:n))**(params%p/2.0))
          end if

       end select

     end subroutine evaltensor_f

     subroutine evaltensor_J(status, n, m, s, J, params)
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(out)   :: J
       class( params_base_type ), intent(in) :: params

       integer :: ii, jj, kk

       ! The function we need to return is 
       !  g_i + H_i s 
       ! where h_i and H_i are the gradient and hessians of the original problem

       select type(params)
       type is(tensor_params_type)
          ! note: params%m contains 'm' from the original problem
          ! if we're passing in the reg. factor via the function/Jacobian, then 
          ! 'm' here is m+n from the original problem
          J(1:n*m) = zero
          ! first, copy in the original J...
!          J(1:n*m) = params%J(1:n*m)
          do jj = 1,n ! columns
             J( (jj-1)*m + 1 : (jj-1)*m + params%m) &
                  = params%J((jj-1)*params%m + 1 : jj*params%m)
             do ii = 1,params%m ! rows 
                do kk = 1,n
                   J( (jj-1)*m + ii) = J( (jj-1)*m + ii) + params%Hi(jj,kk,ii)*s(kk)
                end do
             end do
          end do
          if (params%extra == 1) then 
             ! we're passing in the regularization via the function/Jacobian
             do ii = 1,n ! loop over the columns...
                J(m*(ii-1) + params%m + ii) = sqrt(1.0/params%Delta)
             end do
          elseif (params%extra == 2) then 
             do ii = 1, n
                J(m*ii) = sqrt( (params%p)/(2.0 * params%Delta) ) * & 
                                      (norm2(s(1:n))**( (params%p/2.0) - 2.0)) * & 
                                      s(ii)
             end do
          end if
       end select


     end subroutine evaltensor_J

     subroutine evaltensor_HF(status, n, m, s, f, HF, params)
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(in)   :: f
       real(wp), dimension(*), intent(out) :: HF
       class( params_base_type ), intent(in) :: params

       integer :: ii, jj, kk
       real(wp) :: normx, hf_local
       
       HF(1:n*n) = zero
       select type(params)
       type is (tensor_params_type)
          do ii = 1,params%m
             do jj = 1,n
                do kk = 1,n
                   HF( (jj - 1)*n + kk ) = HF( (jj -1)*n + kk) + & 
                        f(ii)*params%Hi(jj,kk,ii)
                end do
             end do
          end do
          if (params%extra == 2) then 
             normx = norm2(s(1:n))
 !            write(*,*) 'f(m) = ', f(m)
             do jj = 1,n
                do kk = 1,n
                   hf_local = s(jj)*s(kk)
                   if (jj == kk) hf_local = hf_local + normx**2
                   hf_local = sqrt(params%p / (2.0_wp * params%Delta)) * hf_local
                   hf_local = normx**((params%p/2.0_wp) - 4.0_wp) * hf_local
!                   write(*,*) 'hf_local = ',hf_local
                   HF( (jj-1)*n + kk ) = HF( (jj -1)*n + kk) + & 
                        f(m) * hf_local
                end do
             end do
          end if
       end select

       
     end subroutine evaltensor_HF
       

     subroutine get_Hi(n, m, X, params, i, Hi, eval_HF, inform, weights)
       integer, intent(in) :: n, m 
       real(wp), intent(in) :: X(:)
       class( params_base_type ) :: params
       integer, intent(in) :: i 
       real(wp), intent(out) :: Hi(:,:)
       procedure( eval_hf_type ) :: eval_HF
       type( nlls_inform ), intent( inout ) :: inform
       real( wp ), dimension( m ), intent( in ), optional :: weights

       real( wp ) :: ei( m )

       ei = zero
       ei(i) = one

       if ( present(weights) ) then
          call eval_HF(inform%external_return, n, m, X, weights(1:m)*ei, Hi, params)
       else
          call eval_HF(inform%external_return, n, m, X, ei, Hi, params)
       end if
       
       
     end subroutine get_Hi


   end module ral_nlls_internal

