! ral_nlls_double :: a nonlinear least squares solver

module ral_nlls_double

!  use RAL_NLLS_DTRS_double
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

  public :: nlls_solve, nlls_iterate, nlls_finalize
    
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
                         options, inform )
    
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
      
    integer  :: i
    
    type ( NLLS_workspace ) :: w
    
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
       
       call nlls_iterate(n, m, X,                   & 
                         w, &
                         eval_F, eval_J, eval_HF,   & 
                         params,                    &
                         inform, options)
       ! test the returns to see if we've converged
       if (inform%status .ne. 0) then 
          goto 1000 ! error -- exit
       elseif ((inform%convergence_normf == 1).or.(inform%convergence_normg == 1)) then
          goto 1000 ! converged -- exit
       end if
       
     end do main_loop
    
     ! If we reach here, then we're over maxits     
     if (options%print_level > 0 ) write(options%out,1040) 
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
                          inform, options)

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( nlls_inform ), INTENT( OUT ) :: inform
    TYPE( nlls_options ), INTENT( IN ) :: options
    type( NLLS_workspace ), INTENT( INOUT ) :: w
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
      
    integer :: jstatus=0, fstatus=0, hfstatus=0, astatus = 0, svdstatus = 0
    integer :: i, no_reductions, max_tr_decrease = 100
    real(wp) :: rho, normFnew, md, Jmax, JtJdiag
    real(wp) :: FunctionValue, hybrid_tol
    logical :: success, calculate_svd_J
    real(wp) :: s1, sn
    
    ! todo: make max_tr_decrease a control variable

    ! Perform a single iteration of the RAL_NLLS loop
    
    calculate_svd_J = .true. ! todo :: make a control variable 

    if (w%first_call == 1) then
       ! This is the first call...allocate arrays, and get initial 
       ! function evaluations
       if ( options%print_level >= 3 )  write( options%out , 3000 ) 
       ! first, check if n < m
       if (n > m) goto 4070

       call setup_workspaces(w,n,m,options,inform%alloc_status)
       if ( inform%alloc_status > 0) goto 4000

       call eval_F(inform%external_return, n, m, X, w%f, params)
       inform%f_eval = inform%f_eval + 1
       if (inform%external_return .ne. 0) goto 4020

       call eval_J(inform%external_return, n, m, X, w%J, params)
       inform%g_eval = inform%g_eval + 1
       if (inform%external_return .ne. 0) goto 4010
       if (options%relative_tr_radius == 1) then 
          ! first, let's get diag(J^TJ)
          Jmax = 0.0
          do i = 1, n
             ! note:: assumes column-storage of J
             JtJdiag = norm2( w%J( (i-1)*m + 1 : i*m ) )
             if (JtJdiag > Jmax) Jmax = JtJdiag
          end do
          w%Delta = options%initial_radius_scale * (Jmax**2)
          if (options%print_level .ge. 3) write(options%out,3110) w%Delta
       else
          w%Delta = options%initial_radius
       end if

       if ( calculate_svd_J ) then
          call get_svd_J(n,m,w%J,s1,sn,options,svdstatus)
          if ((svdstatus .ne. 0).and.(options%print_level .ge. 3)) then 
             write( options%out, 3000 ) svdstatus
          end if
       end if

       w%normF = norm2(w%f)
       w%normF0 = w%normF

       !    g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g)
       w%g = -w%g
       w%normJF = norm2(w%g)
       w%normJF0 = w%normJF
!       if (options%model == 8 .or. options%model == 9) w%normJFold = normJF
       w%normJFold = w%normJF

       if (options%model == 9) then
          ! make this relative....
          w%hybrid_tol = options%hybrid_tol * ( w%normJF/(0.5*(w%normF**2)) )
       end if
       
       ! save some data 
       inform%obj = 0.5 * ( w%normF**2 )
       inform%norm_g = w%normJF
       inform%scaled_g = w%normJF/w%normF

       if (options%output_progress_vectors) then
          w%resvec(1) = inform%obj
          w%gradvec(1) = inform%norm_g
       end if

       select case (options%model)
       case (1) ! first-order
          w%hf(1:n**2) = zero
       case (2) ! second order
          if ( options%exact_second_derivatives ) then
             call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
             inform%h_eval = inform%h_eval + 1
             if (inform%external_return > 0) goto 4030
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = zero
          end if
       case (9) ! hybrid (MNT)
          ! use first-order method initially
          w%hf(1:n**2) = zero
          w%use_second_derivatives = .false.
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
    if ( calculate_svd_J ) then
       call get_svd_J(n,m,w%J,s1,sn,options,svdstatus)
       if ((svdstatus .ne. 0).and.(options%print_level > 2)) then 
          write( options%out, 3000 ) svdstatus
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

    select case (options%model) ! only update hessians than change..
    case (1) ! first-order
       continue
    case (2) ! second order
       call apply_second_order_info(n,m, & 
            X,w, & 
            eval_Hf, &
            params,options,inform)
       if (inform%external_return .ne. 0) goto 4030
    case (9)
       ! First, check if we need to switch methods
       
       if (w%use_second_derivatives) then 
          if (w%normJF > w%normJFold) then 
             ! switch to Gauss-Newton             
             if (options%print_level .ge. 3) write(options%out,3120) 
             w%use_second_derivatives = .false.
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
             end if
          else 
             w%hybrid_count = 0
          end if
       end if
       
       if (w%use_second_derivatives) then 
          call apply_second_order_info(n,m, &
               X, w, & 
               eval_Hf, &
               params,options,inform)
          if (inform%external_return .ne. 0) goto 4030
       elseif (.not. options%exact_second_derivatives) then 
          ! if exact_second_derivatives are not needed,
          ! call apply_second_order info so that we update the approximation
           call apply_second_order_info(n,m, &
               X,w, &
               eval_Hf, &
               params,options,inform)
          if (inform%external_return .ne. 0) goto 4030
          w%hf(1:n**2) = zero
       else 
          w%hf(1:n**2) = zero
       end if
    end select

    ! update the stats 
    inform%obj = 0.5*(w%normF**2)
    inform%norm_g = w%normJF
    inform%scaled_g = w%normJF/w%normF
    if (options%output_progress_vectors) then
       w%resvec (w%iter + 1) = inform%obj
       w%gradvec(w%iter + 1) = inform%norm_g
    end if
    
    if (options%print_level >=3) write(options%out,3010) inform%obj
    if (options%print_level >=3) write(options%out,3060) w%normJF/w%normF

    !++++++++++++++++++!
    ! Test convergence !
    !++++++++++++++++++!
    call test_convergence(w%normF,w%normJF,w%normF0,w%normJF0,options,inform)
    if (inform%convergence_normf == 1) goto 5000 ! <----converged!!
    if (inform%convergence_normg == 1) goto 5010 ! <----converged!!

    if (options%print_level > 2 ) write(options%out,3100) rho

! Non-executable statements

! print level > 0

!1040 FORMAT(/,'RAL_NLLS failed to converge in the allowed number of iterations')

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
! error returns
4000 continue
    ! generic end of algorithm
        ! all (final) exits should pass through here...
    inform%iter = w%iter
!    if (options%output_progress_vectors) then
    if (allocated(w%resvec)) then
       if(.not. allocated(inform%resvec)) then 
          allocate(inform%resvec(w%iter + 1), stat = astatus)
          if (astatus > 0) call allocation_error(options,'nlls_iterate')
          inform%resvec(1:w%iter + 1) = w%resvec(1:w%iter + 1)
       end if
    end if
    if (allocated(w%gradvec)) then
       if(.not. allocated(inform%gradvec)) then 
          allocate(inform%gradvec(w%iter + 1), stat = astatus)
          if (astatus > 0) call allocation_error(options,'nlls_iterate')
          inform%gradvec(1:w%iter + 1) = w%gradvec(1:w%iter + 1)
       end if
    end if

    return

4010 continue
    ! Error in eval_J
    call eval_error(options,inform,'eval_J')
    goto 4000

4020 continue
    ! Error in eval_J
    call eval_error(options,inform,'eval_F')
    goto 4000

4030 continue
    ! Error in eval_HF
    call eval_error(options,inform,'eval_HF')
    goto 4000

4040 continue 
    ! unsupported choice of model
    if (options%print_level > 0) then
       write(options%error,'(a,i0,a)') 'Error: the choice of options%model = ', &
            options%model, ' is not supported'
    end if
    inform%status = ERROR%UNSUPPORTED_MODEL
   goto 4000

4050 continue 
    ! max tr reductions exceeded
    if (options%print_level > 0) then
       write(options%error,'(a)') 'Error: maximum tr reductions reached'
    end if
    inform%status = ERROR%MAX_TR_REDUCTIONS
    goto 4000

4060 continue 
    ! x makes no progress
    if (options%print_level > 0) then
       write(options%error,'(a)') 'No further progress in X'
    end if
    inform%status = ERROR%X_NO_PROGRESS
    goto 4000

4070 continue
    ! n > m on entry
    if (options%print_level > 0) then
       write(options%error,'(a)') ''
    end if
    inform%status = ERROR%N_GT_M
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

end module ral_nlls_double

