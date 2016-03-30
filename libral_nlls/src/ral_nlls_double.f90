! ral_nlls_double :: a nonlinear least squares solver

module ral_nlls_double

  use ral_nlls_internal

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
    
  public :: nlls_options, nlls_inform, nlls_workspace
  public :: params_base_type
    
contains


  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ S O L V E ****************************!!
  !!******************************************************!!
  !!******************************************************!!

  SUBROUTINE NLLS_SOLVE( n, m, X,                   & 
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
    character (len = 80) :: error_string
    
!!$    write(*,*) 'Controls in:'
!!$    write(*,*) control
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

  subroutine nlls_iterate(n, m, X,                   & 
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
    real(wp) :: rho, normFnew, md, Jmax, JtJdiag
    real(wp) :: FunctionValue
    logical :: success
    
    ! todo: make max_tr_decrease a control variable

    ! Perform a single iteration of the RAL_NLLS loop
    
    if (w%first_call == 1) then
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       !! This is the first call...allocate arrays, and get initial !!
       !! function evaluations                                      !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       if ( options%print_level >= 3 ) write( options%out, 3000 ) 
       if ( options%print_level >= 1 ) write( options%out, 1000 )
       ! first, check if n < m
       if (n > m) goto 4070
       
       ! allocate space for vectors that will be used throughout the algorithm
       call setup_workspaces(w,n,m,options,inform)
       if ( inform%alloc_status > 0) goto 4000

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
               options,svdstatus,w%get_svd_J_ws)
          if ((svdstatus .ne. 0).and.(options%print_level .ge. 3)) then 
             write(options%out,'(a,i0)') 'warning! svdstatus = ', svdstatus
             write( options%out, 3140 ) svdstatus
          end if
       end if

       w%normF = norm2(w%f)
       w%normF0 = w%normF

       !    g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g)
       w%g = -w%g
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
       
       !! Select the order of the model to be used..
       select case (options%model)
       case (1) ! first-order
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
       case (2) ! second order
          if ( options%exact_second_derivatives ) then
             if ( present(weights) ) then
                call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
             else
                call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
             end if
             inform%h_eval = inform%h_eval + 1
             if (inform%external_return > 0) goto 4030
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = zero
          end if
          w%use_second_derivatives = .true.
       case (3) ! hybrid (MNT)
          ! set the tolerance :: make this relative
          w%hybrid_tol = options%hybrid_tol * ( w%normJF/(0.5*(w%normF**2)) )                   
          ! use first-order method initially
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
          if (.not. options%exact_second_derivatives) then 
             ! initialize hf_temp too 
             w%hf_temp(1:n**2) = zero
          end if
       case default
          goto 4040 ! unsupported model -- return to user
       end select
       
       
    end if

    w%iter = w%iter + 1
    if ( options%print_level >= 3 )  write( options%out , 3030 ) w%iter
    inform%iter = w%iter
    
    rho  = -one ! intialize rho as a negative value
    success = .false.
    no_reductions = 0

    do while (.not. success) ! loop until successful
       no_reductions = no_reductions + 1
       if (no_reductions > max_tr_decrease+1) goto 4050

       !+++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the step                        !
       !    d                                      !   
       ! that the model thinks we should take next !
       !+++++++++++++++++++++++++++++++++++++++++++!
       call calculate_step(w%J,w%f,w%hf,w%g,n,m,w%Delta,w%d,w%normd,options,inform,& 
            w%calculate_step_ws)
       if (inform%status .ne. 0) goto 4000
       
       !++++++++++++++++++!
       ! Accept the step? !
       !++++++++++++++++++!
       w%Xnew = X + w%d
       call eval_F(inform%external_return, n, m, w%Xnew, w%fnew, params)
       inform%f_eval = inform%f_eval + 1
       if (inform%external_return .ne. 0) goto 4020
       if ( present(weights) ) then
          ! set f -> Wf
          w%fnew(1:m) = weights(1:m)*w%fnew(1:m)
       end if       
       normFnew = norm2(w%fnew)
       
       !++++++++++++++++++++++++++++!
       ! Get the value of the model !
       !      md :=   m_k(d)        !
       ! evaluated at the new step  !
       !++++++++++++++++++++++++++++!
       call evaluate_model(w%f,w%J,w%hf,w%d,md,m,n,options,w%evaluate_model_ws)
       
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       ! Calculate the quantity                                   ! 
       !   rho = 0.5||f||^2 - 0.5||fnew||^2 =   actual_reduction  !
       !         --------------------------   ------------------- !
       !             m_k(0)  - m_k(d)         predicted_reduction !
       !                                                          !
       ! if model is good, rho should be close to one             !
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       call calculate_rho(w%normF,normFnew,md,rho)
       if (rho > options%eta_successful) success = .true.
       
       !++++++++++++++++++++++!
       ! Update the TR radius !
       !++++++++++++++++++++++!
       call update_trust_region_radius(rho,options,inform,w)
       if (inform%status .ne. 0) goto 4000
       
       if (.not. success) then
          ! finally, check d makes progress
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
            options,svdstatus,w%get_svd_J_ws)
       if ((svdstatus .ne. 0).and.(options%print_level >= 3)) then 
          write( options%out, 3140 ) svdstatus
       end if
    end if
    
    ! g = -J^Tf
    call mult_Jt(w%J,n,m,w%f,w%g)
    w%g = -w%g
    
    w%normJFold = w%normJF
    w%normF = normFnew
    w%normJF = norm2(w%g)
    
    ! setup the vectors needed if second derivatives are not available
    if (.not. options%exact_second_derivatives) then 
       
       w%y       = w%g_old   - w%g
       w%y_sharp = w%g_mixed - w%g

    end if
    
    if (options%model == 3) then
       ! hybrid method -- check if we need second derivatives
       
       if (w%use_second_derivatives) then 
          if (w%normJF > w%normJFold) then 
             ! switch to Gauss-Newton             
             if (options%print_level .ge. 3) write(options%out,3120) 
             w%use_second_derivatives = .false.
             ! save hf as hf_temp
             w%hf_temp(1:n**2) = w%hf(1:n**2)
             w%hf(1:n**2) = zero
          end if
       else
          FunctionValue = 0.5 * (w%normF**2)
          if ( w%normJF/FunctionValue < w%hybrid_tol ) then 
             w%hybrid_count = w%hybrid_count + 1
             if (w%hybrid_count == options%hybrid_switch_its) then
                ! use (Quasi-)Newton
                if (options%print_level .ge. 3) write(options%out,3130) 
                w%use_second_derivatives = .true.
                w%hybrid_count = 0
                ! copy hf from hf_temp
                if (.not. options%exact_second_derivatives) then
                   w%hf(1:n**2) = w%hf_temp(1:n**2)
                end if
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
       write(options%out,1010) w%iter, w%Delta, inform%obj, inform%norm_g, inform%scaled_g
    end if

    !++++++++++++++++++!
    ! Test convergence !
    !++++++++++++++++++!
    call test_convergence(w%normF,w%normJF,w%normF0,w%normJF0,options,inform)
    if (inform%convergence_normf == 1) goto 5000 ! <----converged!!
    if (inform%convergence_normg == 1) goto 5010 ! <----converged!!

    if (options%print_level >= 3 ) write(options%out,3100) rho

! Non-executable statements

! print level > 0

!1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

1000 FORMAT('iter',4x,'Delta',9x,'0.5||f||^2',4x,'||J''f||',7x,'||J''f||/||f||')
1010 FORMAT(   i4, 2x,ES12.4,2x,ES12.4,2x,ES12.4,2x,ES12.4)

! print level > 1



! print level > 2
3000 FORMAT(/,'* Running RAL_NLLS *')
3010 FORMAT('0.5 ||f||^2 = ',ES12.4)
3030 FORMAT('== Starting iteration ',i0,' ==')
3060 FORMAT('||J''f||/||f|| = ',ES12.4)

!3090 FORMAT('Step was unsuccessful -- rho =', ES12.4)
3100 FORMAT('Step was successful -- rho =', ES12.4)
3110 FORMAT('Initial trust region radius taken as ', ES12.4)
3120 FORMAT('** Switching to Gauss-Newton **')
3130 FORMAT('** Switching to (Quasi-)Newton **')
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
       if (inform%alloc_status > 0) goto 4080
       inform%resvec(1:w%iter + 1) = w%resvec(1:w%iter + 1)
    end if
    if (allocated(w%gradvec)) then
       if (allocated(inform%gradvec)) deallocate(inform%gradvec)
       allocate(inform%gradvec(w%iter + 1), stat = inform%alloc_status)
       if (inform%alloc_status > 0) goto 4080
       inform%gradvec(1:w%iter + 1) = w%gradvec(1:w%iter + 1)
    end if
    if (options%calculate_svd_J) then
       if (allocated(inform%smallest_sv) ) deallocate(inform%smallest_sv)
       allocate(inform%smallest_sv(w%iter + 1))
       if (inform%alloc_status > 0) goto 4080
       if (allocated(inform%largest_sv) ) deallocate(inform%largest_sv)
       allocate(inform%largest_sv(w%iter + 1))
       if (inform%alloc_status > 0) goto 4080
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

! convergence 
5000 continue
    ! convegence test satisfied
    if (options%print_level > 2) then
       write(options%out,'(a,i0)') 'RAL_NLLS converged (on ||f|| test) at iteration ', &
            w%iter
    end if
    goto 4000

5010 continue
    if (options%print_level > 2) then
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
    
    w%first_call = 1

    call remove_workspaces(w,options)   
    
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
    else 
       inform%error_message = 'Unknown error number'           
    end if
    
  end subroutine nlls_strerror
      
  
end module ral_nlls_double

