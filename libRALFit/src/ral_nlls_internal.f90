! ral_nlls_internal :: internal subroutines for ral_nlls

module ral_nlls_internal

  use RAL_NLLS_DTRS_double
  use RAL_NLLS_DRQS_double
  use ral_nlls_workspaces
  Use ral_nlls_printing

  implicit none

  private

  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(out) :: f
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
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
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hf_type
  end interface

  abstract interface
     subroutine eval_hp_type(status, n, m, x, y, hp, params)
       import :: params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       double precision, dimension(*), intent(in)  :: x
       double precision, dimension(*), intent(in)  :: y
       double precision, dimension(*), intent(out) :: hp
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hp_type
  end interface

    public :: nlls_solve, nlls_iterate, nlls_finalize, nlls_strerror
    public :: solve_galahad, findbeta, mult_j
    public :: mult_jt, solve_spd, solve_general, matmult_inner
    public :: matmult_outer, outer_product, min_eig_symm, max_eig, all_eig_symm
    public :: remove_workspaces, setup_workspaces
    public :: calculate_step, evaluate_model
    public :: update_trust_region_radius, apply_second_order_info, rank_one_update
    public :: test_convergence, calculate_rho
    public :: solve_LLS, shift_matrix
    public :: dogleg, more_sorensen, generate_scaling, solve_newton_tensor, aint_tr
    public :: switch_to_quasi_newton

contains

  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ S O L V E ****************************!!
  !!******************************************************!!
  !!******************************************************!!

  RECURSIVE SUBROUTINE NLLS_SOLVE( n, m, X,                   &
                         eval_F, eval_J, eval_HF,             &
                         params,                              &
                         options, inform, weights, eval_HP )

!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees,
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments
    implicit none

    INTEGER, INTENT( IN ) :: n, m
    REAL( wp ), DIMENSION( n ), INTENT( INOUT ) :: X
    TYPE( NLLS_inform ), INTENT( OUT ) :: inform
    TYPE( NLLS_options ), INTENT( IN ) :: options
    procedure( eval_f_type ) :: eval_F
    procedure( eval_j_type ) :: eval_J
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
    real( wp ), dimension( m ), intent(in), optional :: weights
    procedure( eval_hp_type ), optional :: eval_HP

    integer  :: i, nrec
    Character(Len=80) :: rec(3)

    type ( NLLS_workspace ) :: w
    Type ( NLLS_workspace ), target :: inner_workspace
!   Link the inner_workspace to the main workspace
    w%iw_ptr => inner_workspace
!   Self reference inner workspace so recursive call does not fail
    inner_workspace%iw_ptr => inner_workspace
    ! Check some user options
    Call Check_options(options, inform)
    If (inform%status/=0) Then
      Go To 100
    End If
    If (buildmsg(1, .False., options)) Then
      Write(rec(1),Fmt=6000)
      Write(rec(2),Fmt=6001)
      Write(rec(3),Fmt=6000)
      nrec = 3
      Call printmsg(1, .False., options, nrec, rec)
    End If

!   Print options, (default no)
    If (buildmsg(0, .False., options).And.options%print_options) Then
      Call print_options(options)
    End If

    main_loop: do i = 1,options%maxit
       if ( present(weights) ) then
          if ( present(eval_HP) ) then
             call nlls_iterate(n, m, X,      &
                  w,                         &
                  eval_F, eval_J, eval_HF,   &
                  params,                    &
                  inform, options, weights=weights,eval_HP=eval_HP)
          else
             call nlls_iterate(n, m, X,      &
                  w,                         &
                  eval_F, eval_J, eval_HF,   &
                  params,                    &
                  inform, options, weights=weights)
          end if
       else
          if ( present(eval_HP) ) then
             call nlls_iterate(n, m, X,      &
                  w,                         &
                  eval_F, eval_J, eval_HF,   &
                  params,                    &
                  inform, options, eval_HP=eval_HP)
          else
             call nlls_iterate(n, m, X,      &
                  w,                         &
                  eval_F, eval_J, eval_HF,   &
                  params,                    &
                  inform, options)
          end if
       end if

       ! test the returns to see if we've converged
       if (inform%status /= 0) then
          ! error -- exit
          goto 100
       elseif ((inform%convergence_normf == 1).or.&
               (inform%convergence_normg == 1).or.&
               (inform%convergence_norms == 1)) then
          ! converged -- exit
          goto 100
       end if
     end do main_loop

     ! If we reach here, then we're over maxits
     inform%status = NLLS_ERROR_MAXITS

100 continue

     call nlls_finalize(w,options)
     If (inform%status /= 0) then
       call nlls_strerror(inform)
     End If
     If (buildmsg(1,.False.,options)) then
       Call print_bye(options,inform)
     End If

6000 Format(1X,57('-'))
6001 Format(2X,'RALFit: An unconstrained nonlinear least-squares solver')
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
                          inform, options, weights, eval_HP)
    implicit none
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
    procedure( eval_hp_type ), optional :: eval_HP

    integer :: i, no_reductions, max_tr_decrease, prncnt
    real(wp) :: rho, rho_gn, normFnew, normJFnew, md, md_gn, Jmax, JtJdiag
    real(wp) :: FunctionValue, normX
    logical :: success
    logical :: bad_allocate
    character :: second
    integer :: num_successful_steps
    integer :: nrec, ierr_dummy
    Character(Len=100) :: rec(3)
    Character(Len=1) :: it_type, inn_flag

!   thread-safe inits
    max_tr_decrease = 100
    bad_allocate = .false.
    num_successful_steps = 0
    ! it_type: Iteration type: R = regular, I = inner
    it_type = 'R'
    Select Type(params)
    Type Is (tensor_params_type)
      it_type = 'I'
    End Select
    ! inn_flag: status flag for the success of inner it convergence
    ! Three-state '-' in inner iteration, 'C' Converged, 'E' not convErged
    inn_flag = '-'
    ! todo: make max_tr_decrease a control variable

    ! Perform a single iteration of the RAL_NLLS loop
    if (w%first_call == 1) then
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       !! This is the first call...allocate arrays, and get initial !!
       !! function evaluations and see if problem is already solved !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
       ! first, check if n < m
       If (n > m) Then
         inform%status = NLLS_ERROR_N_GT_M
         goto 100
       End If
       ! set scalars...
       w%first_call = 0
       w%tr_nu = options%radius_increase
       w%tr_p = 7
       inform%status = 0
       inform%iter = 0
       inform%external_return = 0
       ! allocate space for vectors that will be used throughout the algorithm
       if (options%setup_workspaces) then
          call setup_workspaces(w,n,m,options,inform)
          if ( inform%status/=0) goto 100
       elseif (.not. w%allocated) then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          goto 100
       end if

       ! if needed, setup assign eval_Hp
       if (options%model==4 .and. present(eval_HP)) then
          w%calculate_step_ws%solve_newton_tensor_ws%tparams%eval_HP => eval_HP
          w%calculate_step_ws%solve_newton_tensor_ws%tparams%eval_hp_provided = .true.
       else
          w%calculate_step_ws%solve_newton_tensor_ws%tparams%eval_hp_provided = .false.
       end if

       ! evaluate the residual
       call eval_F(inform%external_return, n, m, X, w%f, params)
       inform%f_eval = inform%f_eval + 1
       If (inform%external_return /= 0) Then
         inform%external_name = 'eval_F'
         inform%status = NLLS_ERROR_EVALUATION
         goto 100
       End If
       if ( present(weights)) then
          ! set f -> Wf
          w%f(1:m) = weights(1:m)*w%f(1:m)
       end if

       ! and evaluate the jacobian
       call eval_J(inform%external_return, n, m, X, w%J, params)
       inform%g_eval = inform%g_eval + 1
       If (inform%external_return /= 0) Then
         inform%external_name = 'eval_J'
         inform%status = NLLS_ERROR_EVALUATION
         goto 100
       End If
       if ( present(weights) ) then
          call scale_J_by_weights(w%J,n,m,weights,options)
       end if

       if (options%relative_tr_radius == 1) then
          ! first, let's get diag(J^TJ)
          Jmax = 0.0_wp
          do i = 1, n
             if (options%Fortran_Jacobian) then
                JtJdiag = norm2( w%J( (i-1)*m + 1 : i*m ) )
             else
                JtJdiag = norm2( w%J( i : n*m : n ) )
             end if
             if (JtJdiag > Jmax) Jmax = JtJdiag
          end do
          w%Delta = options%initial_radius_scale * (Jmax**2)
       else
          w%Delta = options%initial_radius
       end if

       w%normF = norm2(w%f(1:m))
       if (options%regularization > 0 ) then
          normX = norm2(X(1:n))
          call update_regularized_normF(w%normF,normX,options)
       end if
       w%normF0 = w%normF

       !    g = -J^Tf
       call mult_Jt(w%J,n,m,w%f,w%g,options)
       w%g = -w%g
       if (options%regularization > 0) call update_regularized_gradient(w%g,X,normX,options)
       w%normJF = norm2(w%g)
       w%normJF0 = w%normJF
       w%normJFold = w%normJF

       ! save some data
       inform%obj = 0.5_wp * ( w%normF**2 )
       inform%norm_g = w%normJF
       inform%scaled_g = w%normJF/w%normF

       ! if we need to output vectors of the history of the residual
       ! and gradient, the set the initial values
       if (options%output_progress_vectors) then
          w%resvec(1) = inform%obj
          w%gradvec(1) = inform%norm_g
       end if

       !++++++++++++++++++!
       ! Test convergence !
       !++++++++++++++++++!
       ! Note: pretty printing was pushed into test_convergence
       ! norm_2_d is not yet available, test on step size is ignored
       call test_convergence(w%normF,w%normJF,w%normF0,w%normJF0,1.0_wp,options,inform)
       if (inform%convergence_normf == 1 .Or. inform%convergence_normg == 1) Then
!        Converged!
         Go To 100
       End If

       ! set the reg_order to that in the options
       w%calculate_step_ws%reg_order = options%reg_order

       !! Select the order of the model to be used..
       select case (options%model)
       case (1) ! first-order
          w%hf(1:n**2) = 0.0_wp
          w%use_second_derivatives = .false.
          if ( ( options%type_of_method == 2) .and. &
               (options%reg_order .le. 0.0_wp)) then
             ! regularization method, use optimal reg
             w%calculate_step_ws%reg_order = 2.0_wp
          end if
       case (2) ! second order
          if ( options%exact_second_derivatives ) then
             if ( present(weights) ) then
                call eval_HF(inform%external_return, n, m, X, weights(1:m)*w%f, w%hf, params)
             else
                call eval_HF(inform%external_return, n, m, X, w%f, w%hf, params)
             end if
             inform%h_eval = inform%h_eval + 1
             If (inform%external_return /= 0) Then
               inform%external_name = 'eval_HF'
               inform%status = NLLS_ERROR_EVALUATION
               goto 100
             End If

             if (options%regularization > 0) then
                call update_regularized_hessian(w%hf,X,n,options)
             end if
          else
             ! S_0 = 0 (see Dennis, Gay and Welsch)
             w%hf(1:n**2) = 0.0_wp
          end if
          w%use_second_derivatives = .true.
          if ( ( options%type_of_method == 2) .and. &
               (options%reg_order .le. 0.0_wp)) then
             ! regularization method, use optimal reg
             w%calculate_step_ws%reg_order = 3.0_wp
          end if
       case (3) ! hybrid (MNT)
          ! set the tolerance :: make this relative
          w%hybrid_tol = options%hybrid_tol * ( w%normJF/(0.5_wp*(w%normF**2)) )
          ! use first-order method initially
          w%hf(1:n**2) = 0.0_wp
          w%use_second_derivatives = .false.
          if ( (options%type_of_method == 2) .and. &
               (options%reg_order .le. 0.0_wp)) then
             ! regularization method, use optimal reg
             w%calculate_step_ws%reg_order = 2.0_wp
          end if
          if (.not. options%exact_second_derivatives) then
             ! initialize hf_temp too
             w%hf_temp(1:n**2) = 0.0_wp
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
             If (inform%external_return /= 0) Then
               inform%external_name = 'eval_HF'
               inform%status = NLLS_ERROR_EVALUATION
               goto 100
             End If
          else
             ! no second derivatives in tensor model
             inform%status = NLLS_ERROR_NO_SECOND_DERIVATIVES
             Go To 100
          end if
       case default
          inform%status = NLLS_ERROR_UNSUPPORTED_MODEL
          goto 100
       end select

       rho  = -1.0_wp ! intialize rho as a negative value

!      @TAG:ITER0_BANNER
!      The following section builds and prints the initial
!      banner containing the information related to the initial iterate,
!      it builds different banners for each print_level.
!      Details cand be found in tag ITER_BANNER
       If (it_type=='R') Then
         If (buildmsg(2, .True., options)) Then
           Write(rec(1), Fmt=8002)
           Write(rec(2), Fmt=9002)
           Write(rec(3), Fmt=8002)
           Call printmsg(2, .True., options, 3, rec)
         Else If (buildmsg(3, .False., options)) Then
           Write(rec(1), Fmt=8000)
           Write(rec(2), Fmt=9000)
           Write(rec(3), Fmt=8000)
           Call printmsg(3, .False., options, 3, rec)
         End If
       Else
         If (buildmsg(4, .False., options)) Then
           Write(rec(1), Fmt=8000)
           Call printmsg(4, .False., options, 1, rec)
         End If
       End If

       If (buildmsg(2, .True., options).And.it_type == 'R') Then
         Write(rec(1),Fmt=9010) 0, inform%obj, inform%norm_g,      &
           inform%scaled_g
         nrec = 1
         Call printmsg(2, .True., options, nrec, rec)
       ElseIf (buildmsg(3, .False., options)) Then
         If (it_type == 'R') Then
           Write(rec(1),Fmt=9010) 0, inform%obj, inform%norm_g,      &
             inform%scaled_g, w%Delta, rho, '---'//inn_flag
           nrec = 1
           Call printmsg(3, .False., options, nrec, rec)
         Else
           If (buildmsg(4, .False., options)) Then
             Write(rec(1),Fmt=9010) 0, inform%obj, inform%norm_g,      &
               inform%scaled_g, w%Delta, rho, '--'//it_type//inn_flag
             nrec = 1
             Call printmsg(4, .False., options, nrec, rec)
           End If
         End If
       End If

       if (.not. options%exact_second_derivatives) then
          ! let's update g_old, which is needed for call to
          ! rank_one_update
          w%g_old = w%g
       end if

    end if

    w%iter = w%iter + 1
    inform%iter = w%iter
    success = .false.
    no_reductions = 0

    do while (.not. success)
       no_reductions = no_reductions + 1
       If (no_reductions > max_tr_decrease+1) Then
          ! max tr reductions exceeded
          inform%status = NLLS_ERROR_MAX_TR_REDUCTIONS
          goto 100
       End If

       If (buildmsg(3,.False.,options)) Then
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
       call calculate_step(w%J,w%f,w%hf,w%g,&
            X,md,md_gn,&
            n,m,w%use_second_derivatives,w%Delta,eval_HF, params, &
            num_successful_steps, &
            w%Xnew,w%d,w%norm_2_d,w%norm_S_d, &
            options,inform,&
            w%calculate_step_ws, w%tenJ, w%iw_ptr)
       if (inform%status /= 0) then
          if ( (w%use_second_derivatives) .and. &
               (options%model == 3) ) then
             ! don't trust this model -- switch to GN
             ! (See Dennis, Gay and Walsh (1981), Section 5)
             call switch_to_gauss_newton(w,n,options)
             w%hf_temp(:) = 0.0_wp
             cycle
          else
             goto 100
          end if
       end if
       ! Update inn_flag
       If (it_type == 'R'.And.options%model==4) Then
         inn_flag = merge('C', 'E', inform%inner_iter_success)
       End If

       !++++++++++++++++++!
       ! Accept the step? !
       !++++++++++++++++++!
       call eval_F(inform%external_return, n, m, w%Xnew, w%fnew, params)
       inform%f_eval = inform%f_eval + 1

       If (inform%external_return .ne. 0) Then
          inform%external_name = 'eval_F'
          inform%status = NLLS_ERROR_EVALUATION
          goto 100
       End If
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
                w%hf_temp(:) = 0.0_wp
                cycle
             end if
          end if
       else
          ! rho seems to be good -- calculate the Jacobian
          if (.not. options%exact_second_derivatives) then
             ! save the value of g_mixed, which is needed for
             ! call to rank_one_update
             ! g_mixed = -J_k^T r_{k+1}
             call mult_Jt(w%J,n,m,w%fnew,w%g_mixed,options)
             w%g_mixed = -w%g_mixed
          end if

          ! evaluate J and hf at the new point
          call eval_J(inform%external_return, n, m, w%Xnew(1:n), w%J, params)
          inform%g_eval = inform%g_eval + 1
          If (inform%external_return /= 0) Then
            inform%external_name = 'eval_J'
            inform%status = NLLS_ERROR_EVALUATION
            goto 100
          End If
          if ( present(weights) ) then
             ! set J -> WJ
             call scale_J_by_weights(w%J,n,m,weights,options)
          end if

          ! g = -J^Tf
          call mult_Jt(w%J,n,m,w%fnew,w%g,options)
          w%g = -w%g
          if ( options%regularization > 0 ) call update_regularized_gradient(w%g,w%Xnew,normX,options)

          normJFnew = norm2(w%g)

          if ( (log(normJFnew)>100.0_wp) .or. (normJFnew/=normJFnew) ) then
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1),Fmt=3120) normJFnew
               nrec=1
               Call Printmsg(5,.False.,options,nrec,rec)
             End If
             rho = -2.0_wp
             ! reset J
             call eval_J(inform%external_return, n, m, X(1:n), w%J, params)
             inform%g_eval = inform%g_eval + 1
             If (inform%external_return /= 0) Then
               inform%external_name = 'eval_J'
               inform%status = NLLS_ERROR_EVALUATION
               goto 100
             End If
             if ( present(weights) ) then
                ! set J -> WJ
                call scale_J_by_weights(w%J,n,m,weights,options)
             end if
             ! reset g
             if (.not. options%exact_second_derivatives) then
                ! this is already saved...
                w%g = w%g_old
             else
                call mult_Jt(w%J,n,m,w%f,w%g,options)
                w%g = -w%g
             end if
          else
             ! success!!
             w%normJFold = w%normJF
             w%normF = normFnew
             w%normJF = normJFnew

             num_successful_steps = num_successful_steps + 1
             success = .true.
          end if
       end if

       !++++++++++++++++++++++!
       ! Update the TR radius !
       !++++++++++++++++++++++!
       call update_trust_region_radius(rho,options,inform,w)
       if (inform%status /= 0) goto 100

       if (.not. success) then
          If (buildmsg(3, .False., options)) Then
            If (it_type=='R') Then
              Write(rec(1),Fmt=9020) inform%iter, w%Delta, rho,                &
                'U'//second//it_type//inn_flag,inform%inner_iter
              nrec = 1
              Call printmsg(3, .False., options, nrec, rec)
            Else
              If (buildmsg(4, .False., options)) Then
                Write(rec(1),Fmt=9020) inform%iter, w%Delta, rho,                &
                  'U'//second//it_type//inn_flag
                nrec = 1
                Call printmsg(4, .False., options, nrec, rec)
              End If
            End If
          End If
!         ! finally, check d makes progress
!         if ( norm2(w%d) < epsmch * norm2(w%Xnew) ) then
!            write(*,*) 'rhs = ', epsmch * norm2(w%Xnew)
!            inform%obj = 0.5*(w%normF**2)
!            ! x makes no progress
!            inform%status = NLLS_ERROR_X_NO_PROGRESS
!            goto 100
!         end if
       end if
    end do

    ! if we reach here, a successful step has been found
    ! update X and f
    X(1:n) = w%Xnew(1:n)
    w%f(1:m) = w%fnew(1:m)

    if (options%model == 3) then
       ! hybrid method -- check if we need second derivatives

       if (w%use_second_derivatives) then
          if (w%normJF > w%normJFold) then
             ! switch to Gauss-Newton
             call switch_to_gauss_newton(w,n,options)
          end if
       else
          FunctionValue = 0.5_wp * (w%normF**2)
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
             call rank_one_update(w%hf_temp,w,n,options)
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
       If (inform%external_return /= 0) Then
         inform%external_name = 'eval_HF'
         inform%status = NLLS_ERROR_EVALUATION
         goto 100
       End If
    end if

    if (.not. options%exact_second_derivatives) then
       ! now let's update g_old, which is needed for
       ! call to rank_one_update
       w%g_old = w%g
    end if

    ! update the stats
    inform%obj = 0.5_wp*(w%normF**2)
    inform%norm_g = w%normJF
    inform%scaled_g = w%normJF/w%normF
    inform%step = w%norm_2_d
    if (options%output_progress_vectors) then
       w%resvec (w%iter + 1) = inform%obj
       w%gradvec(w%iter + 1) = inform%norm_g
    end if

!   @TAG:ITER_BANNER
!   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!   !       Print iterate one-line information (and possibly header)       !
!   !                                                                      !
!   ! Description of columns                                               !
!   ! * Iter:  current iteration counter valuei                            !
!   ! * error: objective value (F)                                         !
!   ! * grad:  norm-two of gradient (JF)                                   !
!   ! * rel grad: relative gradient norm-two (F/JF)                        !
!   ! * Delta: TR radius                                                   !
!   ! * rho: quantity: actual_reduction / predicted_reduction              !
!   ! * S2IF: flags:                                                       !
!   !     * S=TR iteration was successful (S) or unsuccessful (U)          !
!   !     * 2=iteration used 2nd order information (Y) or not (N)          !
!   !     * I=iteration type: regular (R) or inner (I)                     !
!   !     * F=exit flag from inner solver. Has three states:               !
!   !       Subproblem converged (C), or                                   !
!   !       Subproblem not solved (E) or                                   !
!   !       Currently inside subproblem, in which case (-)                 !
!   !     * Note: dashes (-) in the positions of flags `S`, `2` and 'I'    !
!   !       indicate information is not available.                         !
!   ! * inn it: inner iteration cummulative counter                        !
!   ! * step: size of last step                                            !
!   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    prncnt = inform%inner_iter + w%iter
    If ( options%model == 4 ) then
      ! account for last inner iteration
      prncnt = prncnt + 1
    End If
    If (buildmsg(2, .True., options).And.it_type=='R') Then
      ! Level 2: Print banner every k-th
      If (mod(inform%iter, options%print_header)==0) Then
        Write(rec(1), Fmt=8002)
        Write(rec(2), Fmt=9002)
        Write(rec(3), Fmt=8002)
        Call printmsg(2, .True., options, 3, rec)
      End If
    Else If (buildmsg(5, .False., options)) Then
!       Level 5: Always print banner
        Write(rec(1), Fmt=8000)
        Write(rec(2), Fmt=9000)
        Write(rec(3), Fmt=8000)
        Call printmsg(5, .False., options, 3, rec)
    Else If (buildmsg(4, .False., options)) Then
!     Level 4: Print banner every k-th including inner iteration
      If (mod(prncnt, options%print_header)==0) Then
        Write(rec(1), Fmt=8000)
        Write(rec(2), Fmt=9000)
        Write(rec(3), Fmt=8000)
        Call printmsg(4, .False., options, 3, rec)
      End If
     Else If (buildmsg(3, .False., options).And.it_type=='R') Then
!     Level 3: Same as level 2 but long banner version
      If (mod(inform%iter, options%print_header)==0) Then
        Write(rec(1), Fmt=8000)
        Write(rec(2), Fmt=9000)
        Write(rec(3), Fmt=8000)
        Call printmsg(3, .False., options, 3, rec)
      End If
    End If
    If (buildmsg(2, .True., options).And.it_type == 'R') Then
      Write(rec(1), Fmt=9010) inform%iter, inform%obj, inform%norm_g,        &
        inform%scaled_g
      nrec = 1
      Call printmsg(2, .True., options, nrec, rec)
    ElseIf (buildmsg(3, .False., options)) Then
      If (it_type=='R') Then
        Write(rec(1), Fmt=9010) inform%iter, inform%obj, inform%norm_g,        &
          inform%scaled_g, w%Delta, rho, 'S'//second//it_type//inn_flag,       &
          inform%inner_iter, inform%step
        nrec = 1
        Call printmsg(3, .False., options, nrec, rec)
      Else
        If (buildmsg(4, .False., options)) Then
          Write(rec(1), Fmt=9011) inform%iter, inform%obj, inform%norm_g,        &
            inform%scaled_g, w%Delta, rho, 'S'//second//it_type//inn_flag,       &
            inform%step
          nrec = 1
          Call printmsg(4, .False., options, nrec, rec)
        End If
      End If
    End If

    !++++++++++++++++++!
    ! Test convergence !
    !++++++++++++++++++!
    ! Note: pretty printing was pushed into test_convergence
    call test_convergence(w%normF,w%normJF,w%normF0,w%normJF0,w%norm_2_d,options,inform)
    if (inform%convergence_normf == 1 .Or. inform%convergence_normg == 1       &
          .Or. inform%convergence_norms == 1) Then
!      Converged!
       Go To 100
    End If

100  Continue
!   generic end of algorithm
!   all (final) exits should pass through here...
    if (allocated(w%resvec)) then
       if( allocated(inform%resvec)) deallocate(inform%resvec, stat=ierr_dummy)
       allocate(inform%resvec(w%iter + 1), stat = inform%alloc_status)
       if (inform%alloc_status == 0) Then
         inform%resvec(1:w%iter + 1) = w%resvec(1:w%iter + 1)
       else
         bad_allocate = .true.
       end if
    end if
    if ( (.not.bad_allocate) .And. allocated(w%gradvec)) then
       if (allocated(inform%gradvec)) deallocate(inform%gradvec, stat=ierr_dummy)
       allocate(inform%gradvec(w%iter + 1), stat = inform%alloc_status)
       if (inform%alloc_status == 0) then
         inform%gradvec(1:w%iter + 1) = w%gradvec(1:w%iter + 1)
       else
         bad_allocate = .true.
       end if
    end if

    if (bad_allocate) then
       if (allocated(inform%resvec)) deallocate(inform%resvec ,stat=ierr_dummy)
       if (allocated(inform%gradvec)) deallocate(inform%gradvec ,stat=ierr_dummy)
       if (allocated(inform%smallest_sv)) deallocate(inform%smallest_sv ,stat=ierr_dummy)
       if (allocated(inform%largest_sv)) deallocate(inform%largest_sv ,stat=ierr_dummy)
       inform%status = NLLS_ERROR_ALLOCATION
       inform%bad_alloc = 'nlls_iterate'
    end if

8000 Format(85('-'))
9000 Format(1X,' Iter |  error   |    grad    |  rel grad  |  Delta  |   rho   |S2IF| inn it|  step')
8002 Format(44('-'))
9002 Format(1X,' Iter |  error   |    grad    |  rel grad')
!    Successfull iteration Regular iteration
9010 Format(I7,1X,Es10.4e2,2(1X,Es12.5e2),1X,2(Es9.2e2,1X),A4,1X,I7,1X,Es7.1e2)
!    Successfull iteration Internal iteration
9011 Format(I7,1X,Es10.4e2,2(1X,Es12.5e2),1X,2(Es9.2e2,1X),A4,1X,7X,1X,Es7.1e2)
!    Unsuccessfull iteration partial information
9020 Format(I7,37X,1X,2(Es9.2e2,1X),A4,1X,I7)
! print level > 2
3110 FORMAT('Initial trust region radius taken as ', ES12.4)
3120 FORMAT('||J^Tf|| = ', ES12.4, ': too large. Reducing TR radius.')
  end subroutine nlls_iterate


  !!******************************************************!!
  !!******************************************************!!
  !!***** N L L S _ F I N A L I Z E **********************!!
  !!******************************************************!!
  !!******************************************************!!

  subroutine nlls_finalize(w,options)
    implicit none
    type( nlls_workspace ) :: w
    type( nlls_options ) :: options
    ! reset all the scalars
    w%first_call = 1
    w%iter = 0
    w%calculate_step_ws%reg_order = 0.0
    w%use_second_derivatives = .false.
    w%hybrid_count = 0
    w%hybrid_tol = 1.0
    w%tr_nu = 2.0
    w%tr_p = 3

    if (options%remove_workspaces) call remove_workspaces(w,options)

  end subroutine nlls_finalize

  subroutine nlls_strerror(inform)
    implicit none
    type( nlls_inform ), intent(inout) :: inform

    if ( inform%status == NLLS_ERROR_MAXITS ) then
       inform%error_message = 'Maximum number of iterations reached'
    elseif ( inform%status == NLLS_ERROR_EVALUATION ) then
       write(inform%error_message,Fmt=5004) &
            'Error code from user-supplied subroutine',trim(inform%external_name), &
            'passed error =', inform%external_return
    elseif ( inform%status == NLLS_ERROR_UNSUPPORTED_MODEL ) then
       inform%error_message = 'Unsupported model passed in options'
    elseif ( inform%status == NLLS_ERROR_FROM_EXTERNAL ) then
       write(inform%error_message,Fmt=5004) &
            'The external subroutine',trim(inform%external_name), &
            'passed error =', inform%external_return
    elseif ( inform%status == NLLS_ERROR_UNSUPPORTED_METHOD ) then
       inform%error_message = 'Unsupported nlls_method passed in options'
    elseif ( inform%status == NLLS_ERROR_ALLOCATION ) then
       write(inform%error_message,Fmt=5002) &
            'Bad allocation of memory in', trim(inform%bad_alloc)
    elseif ( inform%status == NLLS_ERROR_MAX_TR_REDUCTIONS ) then
       inform%error_message = 'The trust region was reduced the maximum number of times'
    elseif ( inform%status == NLLS_ERROR_X_NO_PROGRESS ) then
       inform%error_message = 'No progress made in X'
    elseif ( inform%status == NLLS_ERROR_N_GT_M ) then
       inform%error_message = 'The problem is overdetermined'
    elseif ( inform%status == NLLS_ERROR_BAD_TR_STRATEGY ) then
       inform%error_message = 'Unsupported tr_update_stategy passed in options'
    elseif ( inform%status == NLLS_ERROR_FIND_BETA ) then
       inform%error_message = 'Unable to find suitable scalar in findbeta subroutine'
    elseif ( inform%status == NLLS_ERROR_BAD_SCALING ) then
       inform%error_message = 'Unsupported value of scale passed in options'
    elseif ( inform%status == NLLS_ERROR_WORKSPACE_ERROR ) then
       inform%error_message = 'Error accessing pre-allocated workspace'
    elseif ( inform%status == NLLS_ERROR_UNSUPPORTED_TYPE_METHOD ) then
       inform%error_message = 'Unsupported value of type_of_method passed in options'
    elseif ( inform%status == NLLS_ERROR_DOGLEG_MODEL ) then
       inform%error_message = 'Model not supported in dogleg (nlls_method=1)'
    elseif ( inform%status == NLLS_ERROR_AINT_EIG_IMAG ) then
       inform%error_message = 'All eigenvalues are imaginary (nlls_method=2)'
    elseif ( inform%status == NLLS_ERROR_AINT_EIG_ODD ) then
       inform%error_message = 'Odd matrix sent to max_eig subroutine (nlls_method=2)'
    elseif ( inform%status == NLLS_ERROR_MS_MAXITS ) then
       inform%error_message = 'Maximum iterations reached in more_sorensen (nlls_method=3)'
    elseif ( inform%status == NLLS_ERROR_MS_TOO_MANY_SHIFTS ) then
       inform%error_message = 'Too many shifts taken in more_sorensen (nlls_method=3)'
    elseif ( inform%status == NLLS_ERROR_MS_NO_PROGRESS ) then
       inform%error_message = 'No progress being made in more_sorensen (nlls_method=3)'
    elseif ( inform%status == NLLS_ERROR_NO_SECOND_DERIVATIVES ) then
       inform%error_message = 'Exact second derivatives needed for tensor model'
    elseif ( inform%status == NLLS_ERROR_WRONG_INNER_METHOD ) then
       inform%error_message = 'Unsupported value of inner_method passed in options'
    elseif ( inform%status ==  NLLS_ERROR_PRINT_LEVEL) then
       inform%error_message = 'Illegal value of print_level in options'
    elseif ( inform%status ==  NLLS_ERROR_NOT_IMPLEMENTED) then
       inform%error_message = 'Combination of method/regularization options not yet implemented'
    elseif ( inform%status ==  NLLS_ERROR_UNEXPECTED) then
       inform%error_message = 'Unexpected error occured'
    else
       inform%error_message = 'Unknown error number'
    end if
5002  Format (A,1X,A)
5004  Format (A,1X,A,1X,A,1X,I0)
  end subroutine nlls_strerror


! below are the truly internal subroutines...

  RECURSIVE SUBROUTINE calculate_step(J,f,hf,g,X,md,md_gn,n,m,use_second_derivatives, &
                            Delta,eval_HF,params,num_successful_steps,&
                            Xnew,d,norm_2_d,norm_S_d,options,inform,w,tenJ, inner_workspace)
    implicit none
! -------------------------------------------------------
! calculate_step, find the next step in the optimization
! -------------------------------------------------------

    REAL(wp), intent(in) :: J(:), f(:), hf(:)
    REAL(wp), intent(inout) :: g(:)
    REAL(wp), intent(inout) :: Delta
    REAL(wp), intent(in) :: X(:)
    procedure( eval_hf_type ) :: eval_HF
    class( params_base_type ) :: params
    integer, intent(in) :: num_successful_steps
    integer, intent(in)  :: n, m
    logical, intent(in) :: use_second_derivatives
    real(wp), intent(out) :: d(:), Xnew(:)
    real(wp), intent(out) :: md, md_gn
    real(wp), intent(out) :: norm_2_d,norm_S_d
    TYPE( nlls_options ), INTENT( IN ) :: options
    TYPE( nlls_inform ), INTENT( INOUT ) :: inform
    TYPE( calculate_step_work ) :: w
    Type( tenJ_type ), Intent(InOut) :: tenJ
    Type( NLLS_workspace ), Intent(InOut) :: inner_workspace

    real(wp) :: md_bad
    integer :: i, jj
    logical :: scaling_used
    real(wp) :: normx
    Character(Len=85) :: rec(1)


    if (.not. w%allocated) then
      inform%status = NLLS_ERROR_WORKSPACE_ERROR
      goto 100
    end if

    scaling_used = .false.
    d(1:n) = 0.0_wp
    w%scale = 1.0_wp

    norm_2_d = 1.0e10_wp ! set norm_2_d so that it can't be undefined later

    if ( options%model == 4 ) then
       ! tensor model -- call ral_nlls again
       inform%inner_iter_success = .False.
       call solve_newton_tensor(J, f, eval_HF, X, n, m, Delta, num_successful_steps,&
            d, md, params, options, inform, &
            w%solve_newton_tensor_ws, tenJ, inner_workspace)
       If (buildmsg(4, .False., options)) Then
         Write(rec(1), Fmt=8000)
         Call printmsg(4, .False., options, 1, rec)
       End If
       norm_2_d = norm2(d(1:n)) ! ||d||_D
       Xnew = X + d
       call evaluate_model(f,J,hf,X,Xnew,d,md_bad,md_gn,m,n,options,inform,w%evaluate_model_ws)
       If (inform%status/=0) Then
         Go To 100
       End If
    else
       ! compute the hessian used in the model

!       w%scale = 0.0_wp
       w%extra_scale = 0.0_wp

       ! Set A = J^T J
       call matmult_inner(J,n,m,w%A,options)
       ! add any second order information...
       ! so A = J^T J + HF
       ! C <- A + B (aliasing)
       ! call add_matrices3(w%A,hf,n**2,w%A)
       ! A <- A + B
       call add_matrices(w%A,hf,n**2)

       ! and, now, let's add on a reg parameter, if needed
       select case (options%regularization)
       case (1)
          do i = 1, n
             w%A(i,i) = w%A(i,i) + options%regularization_term
          end do
          w%extra_scale = options%regularization_term
       case (2)
          if ( .not. allocated(w%xxt) ) then
            allocate (w%xxt(n,n), stat=inform%alloc_status)
            if (inform%alloc_status/=0) then
              inform%status = NLLS_ERROR_ALLOCATION
              inform%bad_alloc = 'calculate_step'
              GoTo 100
            end if
          end if
          call outer_product(X,n,w%xxt)
          normx = norm2(X(1:n))
          if (normx > epsmch) then
             ! add term from J^TJ
             w%A(1:n,1:n) = w%A(1:n,1:n) + &
                  ( options%regularization_term * options%regularization_power / 2.0_wp ) * &
                  normx**(options%regularization_power - 4.0_wp) * w%xxt
             ! since there's extra terms in the 'real' J, add these to the scaling
             do i = 1, n
                ! add the square of the entries of last row of the 'real' Jacobian
                w%extra_scale(i) = &
                     (options%regularization_term * options%regularization_power / 2.0_wp ) * &
                     (normx**(options%regularization_power-4.0_wp)) * X(i)**2.0_wp
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
             if (inform%status /= 0) Go To 100
             scaling_used = .true.
          end if
       end if

       IF (scaling_used) then
          do i = 1,n
             w%v(i) = w%v(i) / w%scale(i)
             call dscal(n, 1/w%scale(i), w%A(:,i), 1)
          end do
          do jj = 1, n ! swap order for efficiency
             do i = 1, n
                w%A(i,jj) = w%A(i,jj) / w%scale(i)
             end do
          end do

          ! Let's test how well other methods work...
!          do i = 1, n
!             call dscal(n, 1/1.0_wp, w%A(:,i), 1)
!          end do
       end IF


       ! (Gauss-)/(Quasi-)Newton method -- solve as appropriate...

       if ( options%type_of_method == 1) then
          select case (options%nlls_method)
          case (1) ! Powell's dogleg
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3000) 'dogleg'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call dogleg(J,f,hf,g,n,m,Delta,d,norm_S_d,options,inform,w%dogleg_ws)
             if (inform%status /= 0) Go To 100
          case (2) ! The AINT method
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3000) 'AINT_TR'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call AINT_TR(J,w%A,f,X,w%v,hf,n,m,Delta,d,norm_S_d,options,inform,w%AINT_tr_ws)
             if (inform%status /= 0) Go To 100
          case (3) ! More-Sorensen
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3000) 'More-Sorensen'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call more_sorensen(w%A,w%v,n,m,Delta,d,norm_S_d,options,inform,w%more_sorensen_ws)
             if (inform%status /= 0) Go To 100
          case (4) ! Galahad
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3000) 'Galahad DTRS'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call solve_galahad(w%A,w%v,n,m,Delta,num_successful_steps, &
                  d,norm_S_d,w%reg_order,options,inform,w%solve_galahad_ws)
             if (inform%status /= 0) Go To 100
          case default
             inform%status = NLLS_ERROR_UNSUPPORTED_METHOD
             goto 100
          end select ! nlls_method
       elseif (options%type_of_method == 2) then
          select case (options%nlls_method)
          case (3) ! home-rolled regularization solver
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3020) 'RALFit Solver'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call regularization_solver(w%A,w%v,n,m,Delta,num_successful_steps, &
                  d,norm_S_d,w%reg_order,options,inform,w%regularization_solver_ws)
             if (inform%status /= 0) Go To 100
          case(4) ! Galahad
             If (buildmsg(5,.False.,options)) Then
               Write(rec(1), Fmt=3020) 'Galahad DRQS'
               Call printmsg(5,.False.,options,1,rec)
             End If
             call solve_galahad(w%A,w%v,n,m,Delta,num_successful_steps, &
                  d,norm_S_d,w%reg_order,options,inform,w%solve_galahad_ws)
             if (inform%status /= 0) Go To 100
          case default
             inform%status = NLLS_ERROR_UNSUPPORTED_METHOD
             goto 100
          end select ! nlls_method
       else
          inform%status = NLLS_ERROR_UNSUPPORTED_TYPE_METHOD
          goto 100
       end if ! type_of_method

        ! reverse the scaling on the step
        if ( (scaling_used) ) then
           do i = 1, n
              d(i) = d(i) / w%scale(i)
           end do
           ! recalculate ||d||
           norm_2_d = norm2(d)
        else
           norm_2_d = norm_S_d
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
        If (inform%status/=0) Then
          Go To 100
        End If

     end if

100 Continue

     If (buildmsg(5,.False.,options)) Then
       If (inform%status==0) Then
         Write(rec(1), Fmt=3010)
       Else
         Write(rec(1), Fmt=3030)
       End If
       Call printmsg(5,.False.,options,1,rec)
     End If

3000 FORMAT('*** Solving the trust region subproblem using ',A,' ***')
3010 FORMAT('*** Subproblem solution found ***')
3020 FORMAT('*** Solving the regularized subproblem using ',A,' ***')
3030 FORMAT('*** Error or Subproblem solution NOT found ***')
8000 Format(85('-'))

   END SUBROUTINE calculate_step

   subroutine generate_scaling(J,A,n,m,scale,extra_scale,w,options,inform)
    implicit none
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

     If (.not. w%allocated) Then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       goto 100
     End If

     select case (options%scale)
     case (1,2,3)
        do ii = 1,n
           temp = extra_scale(ii)
           what_scale: if (options%scale == 1) then
              ! use the scaling present in gsl:
              ! scale by W, W_ii = ||J(i,:)||_2^2
              !              temp = temp + norm2(J( (ii-1)*m + 1 : ii*m) )**2
              if (options%Fortran_Jacobian) then
                 temp = temp + dot_product(J( (ii-1)*m + 1 : ii*m), J( (ii-1)*m + 1 : ii*m) )
              else
                 temp = temp + dot_product(J(ii:n*m:n),J(ii:n*m:n))
              end if
           elseif ( options%scale == 2) then
              ! scale using the (approximate) hessian
              do jj = 1,n
                 temp = temp + A(ii,jj)**2
              end do
           elseif ( options%scale == 3) then
              ! a less sophisticated version of scale = 1
              do jj = 1,m
                 call get_element_of_matrix(J,m,jj,ii,Jij)
                 temp = temp + Jij**2
              end do
           end if what_scale
           trim_scale: if (temp < options%scale_min) then
              if (options%scale_trim_min) then
                 temp = options%scale_min
              else
                 temp = 1.0_wp
              end if
           elseif (temp > options%scale_max) then
              if (options%scale_trim_max) then
                 temp = options%scale_max
              else
                 temp = 1.0_wp
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
!!$        if (inform%status .ne. 0) goto 100
!!$        do ii = 1,n
!!$           ! plain version...
!!$           w%diag(n + 1 - ii) = w%tempvec(ii)
!!$        end do
!!$        ! todo : require_increase, test for trimming
     case default
        inform%status = NLLS_ERROR_BAD_SCALING
        return
     end select

100 continue
   end subroutine generate_scaling

   subroutine switch_to_gauss_newton(w, n, options)
     Implicit None
     type (nlls_workspace), intent(inout) :: w
     integer, intent(in) :: n
     type (nlls_options), intent(in) :: options
     Character(Len=80) :: rec(1)

     If (buildmsg(5,.False.,options)) Then
       Write(rec(1),Fmt=99999)
99999 FORMAT('*** Switching to Gauss-Newton ***')
       Call Printmsg(5,.False.,options,1,rec)
     End If
     w%use_second_derivatives = .false.
     if ((options%type_of_method == 2) .and. &
          (options%reg_order .le. 0.0_wp)) then
        ! switch to optimal regularization
        w%calculate_step_ws%reg_order = 2.0_wp
     end if
     ! save hf as hf_temp
     w%hf_temp(1:n**2) = w%hf(1:n**2)
     w%hf(1:n**2) = 0.0_wp
   end subroutine switch_to_gauss_newton

   subroutine switch_to_quasi_newton(w, n, options)
     Implicit None
     type (nlls_workspace), intent(inout) :: w
     integer, intent(in) :: n
     type (nlls_options), intent(in) :: options
     Character(Len=80) :: rec(1)

     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=99999)
99999 FORMAT('** Switching to (Quasi-)Newton **')
       Call Printmsg(5,.False.,options,1,rec)
     End If

     w%use_second_derivatives = .true.
     if ((options%type_of_method == 2).and. &
        (options%reg_order .le. 0.0_wp)) then
        ! switch to optimal regularization
        w%calculate_step_ws%reg_order = 3.0_wp
     end if
     w%hybrid_count = 0
     ! copy hf from hf_temp
     if (.not. options%exact_second_derivatives) then
        w%hf(1:n**2) = w%hf_temp(1:n**2)
     end if
   end subroutine switch_to_quasi_newton

   SUBROUTINE dogleg(J,f,hf,g,n,m,Delta,d,normd,options,inform,w)
! -----------------------------------------
! dogleg, implement Powell's dogleg method
! -----------------------------------------
     Implicit None
     REAL(wp), intent(in) :: J(:), hf(:), f(:), g(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     TYPE( dogleg_work ) :: w

     real(wp) :: alpha, beta
     Integer :: nstep
     Character(Len=20) :: steplabs(3) = (/'Gauss-Newton    ', 'Steepest Descent', 'Dogleg          '/)
     Character(Len=80) :: rec(1)

     nstep = 0

     if (.not. w%allocated ) then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       goto 100
     End If

     !     Jg = J * g
     call mult_J(J,n,m,g,w%Jg,options)

     alpha = norm2(g)**2 / norm2( w%Jg )**2

     w%d_sd = alpha * g;

     ! Solve the linear problem...
     select case (options%model)
     case (1)
        ! linear model...
        call solve_LLS(J,f,n,m,w%d_gn,inform,w%solve_LLS_ws)
        if ( inform%status /= 0 ) goto 100
     case default
        inform%status = NLLS_ERROR_DOGLEG_MODEL
        goto 100
     end select

     if (norm2(w%d_gn) <= Delta) then
        ! Gauss-Newton step
        d = w%d_gn
        nstep = 1
     else if (norm2( alpha * w%d_sd ) >= Delta) then
        ! Steepest Descent step
        d = (Delta / norm2(w%d_sd) ) * w%d_sd
        nstep = 2
     else
       ! Dogleg step
        w%d_sd = alpha * w%d_sd
        w%ghat = w%d_gn - w%d_sd
        call findbeta(w%d_sd,w%ghat,Delta,beta,inform)
        if ( inform%status /= 0 ) goto 100
        d = w%d_sd + beta * w%ghat
        nstep = 3
     end if

     normd = norm2(d)

100 continue
     If (buildmsg(5,.false.,options).And.nstep > 0) Then
         Write(rec(1), Fmt=99999) Trim(steplabs(nstep))
         Call Printmsg(5,.False.,options,1,rec)
     End If
99999 Format(A,1X,'step taken')
   END SUBROUTINE dogleg

   SUBROUTINE AINT_tr(J,A,f,X,v,hf,n,m,Delta,d,normd,options,inform,w)
     ! -----------------------------------------
     ! AINT_tr
     ! Solve the trust-region subproblem using
     ! the method of ADACHI, IWATA, NAKATSUKASA and TAKEDA
     ! -----------------------------------------
     Implicit None
     REAL(wp), intent(in) :: J(:), A(:,:), hf(:), f(:), v(:), X(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( AINT_tr_work ) :: w

     integer :: i, size_hard(2)
     real(wp) :: obj_p0, obj_p1, obj_p0_gn, obj_p1_gn
     REAL(wp) :: norm_p0, tau, lam, eta
     Character(Len=80) :: rec(1)

     If ( .not. w%allocated ) Then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       goto 100
     End If

     ! todo..
     ! seems wasteful to have a copy of A and B in M0 and M1
     ! use a pointer?

     tau = 1.0e-4_wp
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
        w%B(i,i) = 1.0_wp
     end do

     select case (options%model)
     case (1)
        call solve_spd(A,-v,w%LtL,w%p0,n,inform)
        if (inform%status /= 0) goto 100
     case default
        call solve_general(A,-v,w%p0,n,inform,w%solve_general_ws)
        if (inform%status /= 0) goto 100
     end select

     call matrix_norm(w%p0,w%B,norm_p0)

     if (norm_p0 < Delta) then
        ! get obj_p0 : the value of the model at p0
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=2000) 'p0'
          Call Printmsg(5,.False.,options,1,rec)
        End If
        call evaluate_model(f,J,hf,X,X,w%p0,obj_p0,obj_p0_gn,m,n, &
             options, inform, w%evaluate_model_ws)
        if (inform%status /= 0) goto 100
     end if

     w%M0(1:n,1:n) = -w%B
     w%M0(n+1:2*n,1:n) = A
     w%M0(1:n,n+1:2*n) = A
     call outer_product(v,n,w%gtg) ! gtg = Jtg * Jtg^T
     w%M0(n+1:2*n,n+1:2*n) = (-1.0_wp / Delta**2) * w%gtg

     w%M1 = 0.0_wp
     w%M1(n+1:2*n,1:n) = -w%B
     w%M1(1:n,n+1:2*n) = -w%B

     call max_eig(w%M0,w%M1,2*n,lam, w%y, w%y_hardcase, options, inform, w%max_eig_ws)
     if ( inform%status /= 0 ) goto 100

     if (norm2(w%y(1:n)) < tau) then
        ! Hard case
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=2010)
          Call Printmsg(5,.False.,options,1,rec)
        End If
        ! overwrite H onto M0, and the outer prod onto M1...
        size_hard = shape(w%y_hardcase)
        call matmult_outer( matmul(w%B,w%y_hardcase), size_hard(2), n, w%M1_small)
        w%M0_small = A(:,:) + lam*w%B(:,:) + w%M1_small
        ! solve Hq + g = 0 for q
        select case (options%model)
        case (1)
           call solve_spd(w%M0_small,-v,w%LtL,w%q,n,inform)
           if (inform%status /= 0) goto 100
        case default
          call solve_general(w%M0_small,-v,w%q,n,inform,w%solve_general_ws)
          if (inform%status /= 0) goto 100
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! (I think..) and inside...fix


        ! find max eta st ||q + eta v(:,1)||_B = Delta
        call findbeta(w%q,w%y_hardcase(:,1),Delta,eta,inform)
        if ( inform%status /= 0 ) goto 100

        !!!!!      ^^TODO^^    !!!!!
        ! currently assumes B = I !!
        !!!!       fixme!!      !!!!

        w%p1(:) = w%q(:) + eta * w%y_hardcase(:,1)

     else
        select case (options%model)
        case (1)
           call solve_spd(A + lam*w%B,-v,w%LtL,w%p1,n,inform)
           if (inform%status /= 0) goto 100
        case default
           call solve_general(A + lam*w%B,-v,w%p1,n,inform,w%solve_general_ws)
           if (inform%status /= 0) goto 100
        end select
        ! note -- a copy of the matrix is taken on entry to the solve routines
        ! and inside...fix
     end if

     ! get obj_p1: the value of the model at p1
     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=2000) 'p1'
       Call Printmsg(5,.False.,options,1,rec)
     End If
     call evaluate_model(f,J,hf,X,X,w%p1,obj_p1,obj_p1_gn,m,n, &
          options,inform,w%evaluate_model_ws)
     if (inform%status /= 0) goto 100
     ! what gives the smallest objective: p0 or p1?
     if (obj_p0 < obj_p1) then
        d = w%p0
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=2030) 'p0'
          Call Printmsg(5,.False.,options,1,rec)
        End If
     else
        d = w%p1
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=2030) 'p1'
          Call Printmsg(5,.False.,options,1,rec)
        End If
     end if

     normd = norm2(d)

100 continue
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
     ! main output :: d, the soln to the TR subproblem
     ! -----------------------------------------
     Implicit None
     REAL(wp), intent(in) :: A(:,:), v(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: nd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ) :: w

     if (options%use_ews_subproblem) then
        call more_sorensen_ew(A,v,n,m,Delta,d,nd,options,inform,w)
        ! inform%status is passed along
     else
        call more_sorensen_noew(A,v,n,m,Delta,d,nd,options,inform,w)
        ! inform%status is passed along
     end if
   end subroutine more_sorensen

   subroutine more_sorensen_ew(A,v,n,m,Delta,d,nd,options,inform,w)
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
     Implicit None
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
     Character(Len=80) :: rec(2)

     ! The code finds
     !  d = arg min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p|| \leq Delta
     !
     ! set A and v for the model being considered here...

     If (.not. w%allocated) Then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       goto 100
     End If

     local_ms_shift = options%more_sorensen_shift

     ! d = -A\v
     ! in this case, A = AplusSigma
     w%AplusSigma(1:n,1:n) = A(1:n,1:n)

     If (options%force_min_eig_symm) Then
       ! Skip solve_spd_nocopy and jump directly to min_eig_symm
       inform%status = 1
     Else
       call solve_spd_nocopy(w%AplusSigma,-v,d,n,inform)
     End If
     if (inform%status == 0) then
        ! A is symmetric positive definite....
        sigma = 0.0_wp
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6000)
          Call Printmsg(5,.False.,options,1,rec)
        End If
     else
        ! reset the error calls -- handled in the code....
        inform%status = 0
        inform%external_return = 0
        inform%external_name = REPEAT( ' ', 80 )
        call min_eig_symm(A,n,sigma,w%y1,options,inform,w%min_eig_symm_ws)
        if (inform%status /= 0) goto 100
        sigma = -(sigma - local_ms_shift)
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6010)
          Call Printmsg(5,.False.,options,1,rec)
        End If
        ! find a shift that makes (A + sigma I) positive definite,
        ! and solve (A + sigma I) (-v) = d
        call check_shift_and_solve(n,A,w%AplusSigma,v,sigma,d,options,inform,w)
        if (inform%status /= 0) goto 100
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6020)
          Call Printmsg(5,.False.,options,1,rec)
        End If
     end if

     nd = norm2(d)

     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=5000)
       Call Printmsg(5,.False.,options,1,rec)
     End If
     ! now, we're not in the trust region initally, so iterate....
     sigma_shift = 0.0_wp
     no_restarts = 0
     ! set 'small' in the context of the algorithm
     epsilon = max( options%more_sorensen_tol * Delta, options%more_sorensen_tiny )

     ! First, check if we're in the t.r. and adjust accordingly
     if (nd .le. Delta) then
        ! we're within the tr radius
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6030)
          Call Printmsg(5,.False.,options,1,rec)
        End If
        if ( abs(sigma) < options%more_sorensen_tiny ) then
           ! we're good....exit
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6040)
             Call Printmsg(5,.False.,options,1,rec)
           End If
           ! initial point was successful
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=5010) 0, nd, sigma, 0.0_wp
             Write(rec(2), Fmt=5040)
             Call Printmsg(5,.False.,options,2,rec)
           End If
           Go To 100
        else if ( abs( nd - Delta ) < epsilon ) then
           ! also good...exit
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6050)
             Call Printmsg(5,.False.,options,1,rec)
           End If
           ! initial point was successful
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=5010) 0, nd, sigma, 0.0_wp
             Write(rec(2), Fmt=5040)
             Call Printmsg(5,.False.,options,2,rec)
           End If
           Go To 100
        end if
        call findbeta(d,w%y1,Delta,alpha,inform)
        if (inform%status /= 0 ) goto 100
        d = d + alpha * w%y1
        nd = norm2(d)
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=5010) 0, nd, sigma, 0.0_wp
          Call Printmsg(5,.False.,options,1,rec)
        End If
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6060)
          Call Printmsg(5,.False.,options,1,rec)
        End If
        ! also good....exit
        ! initial point was successful
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=5040)
          Call Printmsg(5,.False.,options,1,rec)
        End If
        goto 100
     end if

     do i = 1, options%more_sorensen_maxits
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=5010) i-1, nd, sigma, sigma_shift
          Call Printmsg(5,.False.,options,1,rec)
        End If
        if ( abs(nd  - Delta) .le. epsilon) then
           ! we're within the tr radius -- exit
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6035)
             Call Printmsg(5,.False.,options,1,rec)
           End If
           ! initial point was successful
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=5040)
             Call Printmsg(5,.False.,options,1,rec)
           End If
           goto 100
        end if

        w%q = d ! w%q = R'\d
        CALL DTRSM( 'Left', 'Lower', 'No Transpose', 'Non-unit', n, &
             1, 1.0_wp, w%AplusSigma, n, w%q, n ) ! AplusSigma now holds the chol. factors

        nq = norm2(w%q)
        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=6080) nq
          Call Printmsg(5,.False.,options,1,rec)
        End If

        sigma_shift = ( (nd/nq)**2 ) * ( (nd - Delta) / Delta )
        if (abs(sigma_shift) < max(options%more_sorensen_tiny,epsmch) * abs(sigma) ) then
           if (no_restarts < 1) then
              ! find a shift that makes (A + sigma I) positive definite
              call check_shift_and_solve(n,A,w%AplusSigma,v,sigma,d,options,inform,w)
              if (inform%status /= 0) goto 100
              no_restarts = no_restarts + 1
           else
              ! we're not going to make progress...jump out
              inform%status = NLLS_ERROR_MS_NO_PROGRESS
              goto 100
           end if
        else
          ! don't work if last iteration
          If (i == options%more_sorensen_maxits) Then
            exit
          End If
           sigma = sigma + sigma_shift
           call shift_matrix(A,sigma,w%AplusSigma,n)
           call solve_spd_nocopy(w%AplusSigma,-v,d,n,inform)
        end if
        if (inform%status /= 0) goto 100

        nd = norm2(d)
     end do
     ! maxits reached, not converged
     inform%status = NLLS_ERROR_MS_MAXITS
     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=5010)
       Write(rec(2), Fmt=5020)
       Call Printmsg(5,.False.,options,2,rec)
     End If

100  Continue

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
6035 FORMAT('We''re within the trust region radius')
6040 FORMAT('Sigma tiny, so exit')
6050 FORMAT('||d|| = Delta, so exit')
6060 FORMAT('Return d + alpha*y_1')
6070 FORMAT('Converged! ||d|| = Delta at iteration ',i4)
6080 FORMAT('nq = ',ES12.4)
   end subroutine more_sorensen_ew

   subroutine more_sorensen_noew(A,v,n,m,Delta,d,nd,options,inform,w)
     ! -----------------------------------------
     ! more_sorensen
     ! Solve the trust-region subproblem using
     ! the method of More and Sorensen
     !
     ! Using the implementation as in Algorithm 7.3.4
     ! of Trust Region Methods
     !
     ! main output :: d, the soln to the TR subproblem
     ! -----------------------------------------
     Implicit None
     REAL(wp), intent(in) :: A(:,:), v(:), Delta
     integer, intent(in)  :: n, m
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: nd ! ||d||_D
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ) :: w

     real(wp) :: nq
     real(wp) :: sigma, alpha, sigma_shift
     integer :: i, jj, kk
     real(wp) :: sigma_l, sigma_u
     integer :: region ! 1:N, 2:L, 3:G (2-3:F)
     real(wp) :: normF_A, norminf_A
     real(wp) :: dlange, dasum
     real(wp) :: Aii, abs_Aii_sum, lower_sum, upper_sum, min_Aii, max_lower_sum, max_upper_sum, nv
     real(wp) :: uHu, dHd, theta, kappa_easy, kappa_hard, indef_delta
     integer :: location_of_breakdown
     character(1) :: region_char
     logical :: factorization_done
     Character(Len=80) :: rec(1)

     ! The code finds
     !  d = arg min_p   v^T p + 0.5 * p^T A p
     !       s.t. ||p|| \leq Delta
     !
     ! set A and v for the model being considered here...

     If (.not. w%allocated) Then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       Go To 100
     End If

     theta = 0.5_wp
     kappa_easy = 0.001_wp
     kappa_hard = 0.002_wp

     ! Set up the initial variables....
     normF_A = dlange('F',n,n,A,n,w%norm_work)
     norminf_A = dlange('I',n,n,A,n,w%norm_work)
     do i = 1, n
        Aii = A(i,i)
        abs_Aii_sum = dasum(n,A(:,i),1) ! = sum |A_ij|
        if (Aii < 0.0_wp) then
           lower_sum = abs_Aii_sum + 2*Aii ! A_ij + sum_(i ne j) |A_ij|
           upper_sum = abs_Aii_sum
        else
           lower_sum = abs_Aii_sum
           upper_sum = abs_Aii_sum + 2*Aii ! A_ij + sum_(i ne j) |A_ij|
        end if
        if (i == 1) then
           min_Aii = Aii
           max_lower_sum = lower_sum
           max_upper_sum = upper_sum
        else
           min_Aii = min(Aii,min_Aii)
           max_lower_sum = max(lower_sum, max_lower_sum)
           max_upper_sum = max(upper_sum, max_upper_sum)
        end if
     end do
     nv = norm2(v)

     sigma_l = max(0.0_wp, -min_Aii, nv/Delta - min(max_lower_sum,normF_A,norminf_A))
     sigma_u = max(0.0_wp,nv/Delta + min(max_upper_sum,normF_A,norminf_A))

     sigma = max(0.0_wp,sigma_l)
     ! todo: if the model hasn't changed, use the terminating sigma there instead (p192, TR book)...

     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=6030) sigma_l, sigma, sigma_u
       Call printmsg(5,.False.,options,1,rec)
     End If

     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=5000)
       Call printmsg(5,.False.,options,1,rec)
     End If

     factorization_done = .false.

     ! now start the iteration...
     do i = 1, options%more_sorensen_maxits
        ! d = -A\v
        if (.not. factorization_done) then
           call shift_matrix(A,sigma,w%AplusSigma,n)
           call solve_spd(w%AplusSigma,-v,w%LtL,d,n,inform)
        else
           factorization_done = .false.
        end if
        if (inform%status == 0) then
           ! A is symmetric positive definite....
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6000)
             Call printmsg(5,.False.,options,1,rec)
           End If
           nd = norm2(d)
           if (nd < Delta) then
              region = 3
              region_char = 'G'
              sigma_u = sigma
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6010) sigma_u
                Call printmsg(5,.False.,options,1,rec)
              End If

              ! check for interior convergence....
              if ( abs(sigma) < 1.0e-16_wp) then
                 If (buildmsg(5,.False.,options)) Then
                   Write(rec(1), Fmt=5060)
                   Call printmsg(5,.False.,options,1,rec)
                 End If
                 goto 100
              end if
           else
              region = 2
              region_char = 'L'
              sigma_l = sigma
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6020) sigma_l
                Call printmsg(5,.False.,options,1,rec)
              End If
           end if
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6000)
             Call printmsg(5,.False.,options,1,rec)
           End If
        else
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6050)
             Call printmsg(5,.False.,options,1,rec)
           End If
           location_of_breakdown = inform%external_return
           ! reset the error calls -- handled in the code....
           inform%status = 0
           inform%external_return = 0
           inform%external_name = REPEAT( ' ', 80 )
           region = 1 ! N
           region_char = 'N'
           sigma_l = sigma
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6020) sigma_l
             Call printmsg(5,.False.,options,1,rec)
           End If
        end if

        if (region .ne. 1 ) then ! in F
           w%q = d ! w%q = R'\d
           CALL DTRSM( 'Left', 'Lower', 'No Transpose', 'Non-unit', n, &
                1, 1.0_wp, w%LtL, n, w%q, n )
           nq = norm2(w%q)
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6080) nq
             Call printmsg(5,.False.,options,1,rec)
           End If
           sigma_shift = ( (nd/nq)**2 ) * ( (nd - Delta) / Delta )
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6070) sigma_shift
             Call printmsg(5,.False.,options,1,rec)
           End If
           if (region == 3) then ! lambda in G
              ! use the linpack method...
              call linpack_method(n,w%AplusSigma,w%LtL,w%y1,uHu)
              sigma_l = max(sigma_l, sigma - uHu)
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6020) sigma_l
                Call printmsg(5,.False.,options,1,rec)
              End If
              dHd = dot_product(d, matmul(w%AplusSigma,d))
              call findbeta(d,w%y1,Delta,alpha,inform) ! check -- is this what I need?!?!
              if (inform%status /= 0 ) goto 100
              d = d + alpha * w%y1
              ! check for termination
              if ( (alpha**2) * uHu .le. kappa_hard *(dHd + sigma * Delta**2 )) then
                 If (buildmsg(5,.False.,options)) Then
                   Write(rec(1), Fmt=5050)
                   Call printmsg(5,.False.,options,1,rec)
                 End If
                 goto 100
              end if
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6060)
                Call printmsg(5,.False.,options,1,rec)
              End If
           end if
        else  ! not in F
           ! find an alpha and y1 such that
           !  (A + sigmaI + alpha e_k e_k^T)v = 0
!!$           write(*,*) 'location of breakdown = ', location_of_breakdown
!!$           write(*,*) 'A = '
!!$           do jj = 1, n
!!$              write(*,*) w%AplusSigma(:,jj)
!!$           end do
!!$           write(*,*) 'LtL = '
!!$           do jj = 1, n
!!$              write(*,*) w%LtL(:,jj)
!!$           end do
           indef_delta = -w%AplusSigma(location_of_breakdown,location_of_breakdown)
           do jj = 1, location_of_breakdown-1
              indef_delta = indef_delta + w%LtL(location_of_breakdown,jj)**2
           end do

           w%y1 = 0.0_wp
           w%y1(location_of_breakdown) = 1.0_wp
           do jj = location_of_breakdown - 1, 1, -1
              do kk = jj + 1, location_of_breakdown
                 w%y1(jj) = w%LtL(kk,jj) * w%y1(kk)
              end do
              w%y1(jj) = -w%y1(jj) / w%LtL(jj,jj)
           end do

           sigma_l = max(sigma_l, sigma + indef_delta/norm2(w%y1))
        end if

        ! check for termination
        if ( (region == 2) .or. (region == 3) ) then ! region F
           if (abs(nd - Delta) .le. kappa_easy * Delta) then
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=5030)
                Call printmsg(5,.False.,options,1,rec)
              End If
              goto 100
           end if
           if ( region == 3 ) then ! region G
              if (sigma == 0) then
                 If (buildmsg(5,.False.,options)) Then
                   Write(rec(1), Fmt=5040)
                   Call printmsg(5,.False.,options,1,rec)
                 End If
                 goto 100
              else
!!$                 dHd = dot_product(d, matmul(w%AplusSigma,d))
!!$                 if ( (alpha**2) * uHu .le. kappa_hard *(dHd + sigma * Delta**2 )) then
!!$                    if (options%print_level >= 2) write(options%out,5050)
!!$                    d = d + alpha * w%y1
!!$                    goto 100
!!$                 end if
!!$                 ! update d (as this wasn't done earlier)
!!$                 d = d + alpha * w%y1
              end if
           end if
        end if


        if ( ( region == 2 ) .and. ( norm2(v) > 0) ) then
           sigma = sigma + sigma_shift
            If (buildmsg(5,.False.,options)) Then
              Write(rec(1), Fmt=6040) sigma
              Call printmsg(5,.False.,options,1,rec)
            End If
        elseif (region == 3) then
           call shift_matrix(A,sigma + sigma_shift,w%AplusSigma,n)
           call solve_spd(w%AplusSigma,-v,w%LtL,d,n,inform)
           if (inform%status == 0) then
              region = 2
              sigma = sigma + sigma_shift
              factorization_done = .true.
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6040) sigma
                Call printmsg(5,.False.,options,1,rec)
              End If
           else
              sigma_l = max(sigma_l, sigma+sigma_shift)
              ! check sigma_l for interior convergence
              sigma = max( sqrt(sigma_l * sigma_u), &
                   sigma_l + theta * ( sigma_u - sigma_l) )
              If (buildmsg(5,.False.,options)) Then
                Write(rec(1), Fmt=6090) sigma
                Call printmsg(5,.False.,options,1,rec)
              End If
           end if
        else
           sigma = max( sqrt(sigma_l * sigma_u), &
                sigma_l + theta * ( sigma_u - sigma_l) )
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=6090) sigma
             Call printmsg(5,.False.,options,1,rec)
           End If
        end if

        If (buildmsg(5,.False.,options)) Then
          Write(rec(1), Fmt=5010) i, region_char, nd, sigma_l, sigma, sigma_u
          Call printmsg(5,.False.,options,1,rec)
        End If

     end do

     ! maxits reached, not converged
     inform%status = NLLS_ERROR_MS_MAXITS
     If (buildmsg(5,.False.,options)) Then
       Write(rec(1), Fmt=5020)
       Call printmsg(5,.False.,options,1,rec)
     End If

100 continue

  If (inform%status==0) Then
    nd = norm2(d)
    If (buildmsg(5,.False.,options)) Then
      Write(rec(1), Fmt=5010) i, region_char, nd, sigma_l, sigma, sigma_u
      Call printmsg(5,.False.,options,1,rec)
    End If
  End If

! print_level >= 2
5000 FORMAT('iter',2x,'region',2x,'||d||',9x,'sigma_l',9x,'sigma',9x,'sigma_u')
5010 FORMAT(i4,2x,A,5x,ES12.4,2x,ES12.4,4x,ES12.4,2x,ES12.4)
5020 FORMAT('More-Sorensen failed to converge within max number of iterations')
5030 FORMAT('Converged! | ||d|| - Delta | <= kappa_easy * Delta')
5040 FORMAT('Converged! sigma = 0 in region G')
5050 FORMAT('Converged! MINPACK test')
5060 FORMAT('Converged!')
! print_level >= 3
6000 FORMAT('A + sigma I is symmetric positive definite')
6010 FORMAT('Replacing sigma_u by ',ES12.4)
6020 FORMAT('Replacing sigma_l by ',ES12.4)
6030 FORMAT('sigma_l = ',ES12.4,' sigma = ',ES12.4, ' sigma_u = ', ES12.4)
6040 FORMAT('Replacing sigma by sigma + sigma_shift = ',ES12.4)
6050 FORMAT('A + sigma I is indefinite')
6060 FORMAT('Return d + alpha*y_1')
6070 FORMAT('sigma_shift = ',ES12.4)
6080 FORMAT('nq = ',ES12.4)
6090 FORMAT('Replacing sigma by max(gm,av) = ', ES12.4)
   end subroutine more_sorensen_noew


   subroutine linpack_method(n,H,LLt,u,uHu)

     !-------------------------------------------------
     ! linpack_method
     !
     ! Find a u such that (u,H(sigma)u) is small
     ! Uses the 'linpack method' described in Section 7.3
     ! (page 191) of the Trust Region book
     !-------------------------------------------------

     Implicit None
     integer, intent(in) :: n
     real(wp), intent(in) :: H(:,:), LLt(:,:)
     real(wp), intent(out) :: u(:), uHu

     real(wp) :: v_minus, v_plus, ltimesw
     integer :: i

     u = 0.0_wp

!!$     write(*,*) 'LLt = '
!!$     do i = 1, n
!!$        write(*,*) LLt(:,i)
!!$     end do

     do i = 1, n
        v_plus = 1
        v_minus = -1
        if (i > 1) then
           ltimesw = dot_product(LLt(i,1:i-1),u(1:i-1))
           v_plus = 1 - ltimesw
           v_minus = -1 - ltimesw
        end if
        v_plus = v_plus / LLt(i,i)
        v_minus = v_minus / LLt(i,i)
        if ( v_plus > v_minus) then
           u(i) = v_plus
        else
           u(i) = v_minus
        end if
     end do

     ! in the above, u is the w in the TR book
     ! now set u = L'\u
     CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', n, &
          1, 1.0_wp, LLt, n, u, n )
     u = u / norm2(u)   ! and normalize

!!$     write(*,*) 'H = '
!!$     do i = 1, n
!!$        write(*,*) H(:,i)
!!$     end do

     uHu = dot_product(u, matmul(H,u))
!     write(*,*) 'uHu = ', uHu

   end subroutine linpack_method

   subroutine check_shift_and_solve(n,A,AplusSigma,v,sigma,d,options,inform,w)

     !--------------------------------------------------
     ! check_shift_and_solve
     !
     ! Given an indefinite matrix A, find a ensures that
     ! the input shift sigma makes (A + sigma I) is positive
     ! definite, factorizes (A+sigmaI), and solves (A + sigma I)(-v) = d
     !--------------------------------------------------

     Implicit None
     integer, intent(in) :: n
     real(wp), intent(in) :: A(:,:), v(:)
     real(wp), intent(inout) :: AplusSigma(:,:)
     real(wp), intent(inout) :: sigma
     real(wp), intent(inout) :: d(:)
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform
     type( more_sorensen_work ), intent(inout) :: w

     integer :: no_shifts
     logical :: successful_shift
     Character(Len=80) :: rec(1)

     no_shifts = 0
     successful_shift = .false.
     do while( .not. successful_shift )
        call shift_matrix(A,sigma,w%AplusSigma,n)
        call solve_spd_nocopy(w%AplusSigma,-v,d,n,inform)
        if ( inform%status /= 0 ) then
           ! reset the error calls -- handled in the code....
           inform%status = 0
           inform%external_return = 0
           inform%external_name = REPEAT( ' ', 80 )
           no_shifts = no_shifts + 1
           If ( no_shifts >= 10 ) Then
             inform%status = NLLS_ERROR_MS_TOO_MANY_SHIFTS
             goto 100 ! too many shifts -- exit
           End If
           sigma =  sigma + (10.0_wp**no_shifts) * options%more_sorensen_shift
           If (buildmsg(5,.False.,options)) Then
             Write(rec(1), Fmt=99999) sigma
             Call Printmsg(5,.False.,options,1,rec)
           End If
        else
           successful_shift = .true.
        end if
     end do

100   continue
99999 FORMAT('Trying a shift of sigma = ',ES12.4)
   end subroutine check_shift_and_solve


   subroutine solve_galahad(A,v,n,m,Delta,num_successful_steps,d,normd,reg_order,options,inform,w)

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

     Implicit None
     REAL(wp), intent(in) :: A(:,:), v(:)
     REAL(wp), intent(inout) :: Delta
     integer, intent(in)  :: n, m
     integer, intent(in) :: num_successful_steps
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D, where D is the scaling
     real(wp), intent(in) :: reg_order
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

     If (.not. w%allocated ) Then
       inform%status = NLLS_ERROR_WORKSPACE_ERROR
       goto 100
     End If

     ! We have the unprocessed matrices, we need to get an
     ! eigendecomposition to make A diagonal
     !
     call all_eig_symm(A,n,w%ew,w%ev,w%all_eig_symm_ws,inform)
     if (inform%status /= 0) goto 100

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
           w%v_trans(ii) = 0.0_wp
        end if
        if (abs(w%ew(ii)) < epsmch) then
           w%ew(ii) = 0.0_wp
        end if
     end do

    select case (options%type_of_method)
     case (1)
        call dtrs_initialize( dtrs_options, dtrs_inform )
        ! Does not fail.
        call dtrs_solve(n, Delta, 0.0_wp, w%v_trans, w%ew, w%d_trans, &
                        dtrs_options, dtrs_inform )
        if ( dtrs_inform%status /= 0) then
           inform%external_return = dtrs_inform%status
           inform%external_name = 'galahad_dtrs'
           inform%status = NLLS_ERROR_FROM_EXTERNAL
           goto 100
        end if
     case(2)
        call drqs_initialize( drqs_options, drqs_inform )
        ! Does not fail.

        proceed = .false.
        do while (.not. proceed)
           reg_param = options%base_regularization + 1.0_wp/Delta
           call drqs_solve &
                (n, reg_order, reg_param, 0.0_wp, w%v_trans, w%ew, w%d_trans, &
                drqs_options, drqs_inform )
           if ( drqs_inform%status == -7 ) then
              ! drqs_solve has failed because the matrix
              !     J'J + *1/(2*Delta) * I
              ! is not spd. Fix this by decreasing Delta until it's big enough
              Delta =  Delta / 10.0_wp
              ! alternatively, we could use a higher order regularization:
              ! call drqs_solve(n,3.0_wp,1.0_wp/Delta, 0.0_wp, w%v_trans,&
              !                 w%ew, w%d_trans,drqs_options, drqs_inform)
              continue
           elseif ( drqs_inform%status /= 0) then
              inform%external_return = drqs_inform%status
              inform%external_name = 'galahad_drqs'
              inform%status = NLLS_ERROR_FROM_EXTERNAL
              goto 100
           else
              proceed = .true.
           end if
        end do
     end select
  ! and return the un-transformed vector
     call mult_J(w%ev,n,n,w%d_trans,d)

     normd = norm2(d) ! ||d||_D

100 continue

2000 FORMAT('Regularization order used = ',ES12.4)

   end subroutine solve_galahad

   subroutine regularization_solver(A,v,n,m,Delta,num_successful_steps,&
        d,normd,reg_order,options,inform,w)

     !---------------------------------------------
     ! regularization_solver
     ! Solve the regularized subproblem using
     ! a home-rolled algorithm
     ! Note: inform%status is set by solve_spd
     !--------------------------------------------

     Implicit None
     REAL(wp), intent(in) :: A(:,:), v(:)
     REAL(wp), intent(inout) :: Delta
     integer, intent(in)  :: n, m
     integer, intent(in) :: num_successful_steps
     real(wp), intent(out) :: d(:)
     real(wp), intent(out) :: normd ! ||d||_D, where D is the scaling
     real(wp), intent(in) :: reg_order
     type( regularization_solver_work ) :: w
     TYPE( nlls_options ), INTENT( IN ) :: options
     TYPE( nlls_inform ), INTENT( INOUT ) :: inform

     real(wp) :: reg_param

     reg_param = 1.0_wp/Delta

     if ( reg_order == 2.0_wp ) then
        call shift_matrix(A,reg_param,w%AplusSigma,n)
        call solve_spd(w%AplusSigma,-v,w%LtL,d,n,inform)
        ! informa%status is passed along, this routine exits here
     else
!      Feature not yet implemented, this should have been caught in
!      check_options()
       inform%status = NLLS_ERROR_UNEXPECTED
     end if

     If (inform%status==0) Then
       normd = norm2(d)
     Else
       normd = 1.0e10_wp
     End If
   end subroutine regularization_solver

   SUBROUTINE solve_LLS(J,f,n,m,d_gn,inform,w)

!  -----------------------------------------------------------------
!  solve_LLS, a subroutine to solve a linear least squares problem
!  -----------------------------------------------------------------

       Implicit None
       REAL(wp), DIMENSION(:), INTENT(IN) :: J
       REAL(wp), DIMENSION(:), INTENT(IN) :: f
       INTEGER, INTENT(IN) :: n, m
       REAL(wp), DIMENSION(:), INTENT(OUT) :: d_gn
       type(NLLS_inform), INTENT(INOUT) :: inform

       character(1) :: trans = 'N'
       integer :: nrhs = 1, lwork, lda, ldb
       type( solve_LLS_work ) :: w

       If (.not. w%allocated) Then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       End If

       lda = m
       ldb = max(m,n)
       w%temp(1:m) = f(1:m)
       lwork = size(w%work)

       w%Jlls(:) = J(:)

       call dgels(trans, m, n, nrhs, w%Jlls, lda, w%temp, ldb, w%work, lwork, &
            inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dgels'
          Go To 100
       end if

       d_gn = -w%temp(1:n)

100   continue
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

     Implicit None
     real(wp), dimension(:), intent(in) :: a, b
     real(wp), intent(in) ::  Delta
     real(wp), intent(out) :: beta
     type( nlls_inform ), intent(inout) :: inform

     real(wp) :: c, normb2, norma2, discrim

     c = dot_product(a,b)

     norma2 = norm2(a)**2
     normb2 = norm2(b)**2

     discrim = c**2 + (normb2)*(Delta**2 - norma2);
     if ( discrim < 0.0_wp ) then
        inform%status = NLLS_ERROR_FIND_BETA
        inform%external_name = 'findbeta'
        Go To 100
     end if

     if (c .le. 0.0_wp) then
        beta = (-c + sqrt(discrim) ) / normb2
     else
        beta = (Delta**2 - norma2) / ( c + sqrt(discrim) )
     end if

100  Continue
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

       Implicit None
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

       real(wp) :: xtd, normx, p, sigma
       Character(Len=80) :: rec(1)

       md = 0.0_wp
       md_gn = 0.0_wp

       If (.not. w%allocated ) Then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       End If

       !Jd = J*d
       call mult_J(J,n,m,d,w%Jd,options)

       md_gn = 0.5_wp * norm2(f(1:m) + w%Jd(1:m))**2

       ! if we are solving a regularized problem, update terms
       p = options%regularization_power
       sigma = options%regularization_term
       select case (options%regularization)
       case (1)
          md_gn = md_gn + 0.5_wp * sigma * norm2(Xnew(1:n))**2
       case (2)
          normx = norm2(X(1:n))
          xtd = dot_product(X(1:n),d(1:n))
          md_gn = md_gn + &
               sigma * ( 1.0_wp/p * (normx**p) + &
               (normx**(p-2.0_wp)) * xtd + &
               (p/4.0_wp) * (normx**(p-4.0_wp)) * (xtd**2) )
       end select

       select case (options%model)
       case (1) ! first-order (no Hessian)
          md = md_gn
       case (4) ! tensor model
          ! nothing to do here...
       case default
          ! these have a dynamic H -- recalculate
          ! H = J^T J + HF, HF is (an approx?) to the Hessian
          call mult_J(hf,n,n,d,w%Hd)
          md = md_gn + 0.5_wp * dot_product(d(1:n),w%Hd(1:n))
          ! regularized newton terms taken care of already in apply_second_order_info
       end select
       If (buildmsg(5,.False.,options)) Then
         Write(rec(1), Fmt=99999) md
         Call Printmsg(5,.False.,options,1,rec)
       End If

100   continue
99999 FORMAT('Model evaluated successfully: m_k(d) = ',ES12.4)
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

       Implicit None
       real(wp), intent(in)  :: normf    ! ||f(x_k)|| at
       real(wp), intent(in)  :: normfnew ! ||f(x_k + d)||
       real(wp), intent(in)  :: md       !    m_k(d)
       real(wp), intent(out) :: rho      ! act_red / pred_red (close to 1 == good)
       TYPE( nlls_options ), INTENT( IN ) :: options

       real(wp) :: actual_reduction, predicted_reduction
       real(wp) :: tol
       Character(Len=80) :: rec(3)

       actual_reduction = ( 0.5_wp * (normf**2) ) - ( 0.5_wp * (normfnew**2) )
       predicted_reduction = ( ( 0.5_wp * (normf**2) ) - md )

       tol = 10.0_wp * epsmch

       if ( abs(actual_reduction) < tol ) then
          rho = 1.0_wp
       else if (abs( predicted_reduction ) < tol ) then
          rho = 1.0_wp
       else
          rho = actual_reduction / predicted_reduction
       end if

       If (buildmsg(5,.False.,options)) Then
         Write(rec(1), Fmt=99999) actual_reduction
         Write(rec(2), Fmt=99998) predicted_reduction
         Write(rec(3), Fmt=99997) rho
         Call Printmsg(5,.False.,options,3,rec)
       End If

99999  FORMAT('Actual reduction (in cost function) = ', ES12.4)
99998  FORMAT('Predicted reduction (in model) = ', ES12.4)
99997  FORMAT('rho returned = ', ES12.4)

     end subroutine calculate_rho

     subroutine apply_second_order_info(n,m,X,w,eval_Hf,params,options,inform, &
         weights)
       Implicit None
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
          call rank_one_update(w%hf,w,n,options)
       end if

       ! update the hessian here if we're solving a regularized problem
       if (options%regularization > 0) then
          call update_regularized_hessian(w%hf,X,n,options)
       end if
     end subroutine apply_second_order_info

     subroutine update_regularized_normF(normF,normX,options)
       Implicit None
       real(wp), intent(inout) :: normF
       real(wp), intent(in) :: normX
       type( nlls_options ), intent(in) :: options

       select case(options%regularization)
       case (1)
          normF = sqrt(normF**2 + &
               options%regularization_term * normX**2 )
       case (2)
          normF = sqrt(normF**2 + &
               ( 2.0_wp * options%regularization_term / options%regularization_power )  *  &
               normX**options%regularization_power )
       Case Default
         ! No regularization applied
       end select
     end subroutine update_regularized_normF

     subroutine update_regularized_gradient(g,X,normX,options)
       Implicit None
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
       Implicit None
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
                hf_local = sigma * normx**(p - 4.0_wp) * hf_local
                hf( (ii-1)*n + jj) = hf( (ii-1)*n + jj) + hf_local
             end do
          end do
       end if
     end subroutine update_regularized_hessian

     subroutine rank_one_update(hf,w,n,options)
       Implicit None
       real(wp), intent(inout) :: hf(:)
       type( NLLS_workspace ), intent(inout) :: w
       integer, intent(in) :: n
       type( NLLS_options ), intent(in) :: options

       real(wp) :: yts, alpha, dSks

       w%y       = w%g_old   - w%g
       w%y_sharp = w%g_mixed - w%g

       yts = dot_product(w%d,w%y)
       if ( abs(yts) < 10.0_wp * epsmch ) then
          ! safeguard: skip this update
          Go To 100
       end if

       call mult_J(hf,n,n,w%d,w%Sks) ! hfs = S_k * d

       w%ysharpSks = w%y_sharp - w%Sks

       ! now, let's scale hd (Nocedal and Wright, Section 10.2)
       dSks = abs(dot_product(w%d,w%Sks))
       if ( abs(dSks) < 10.0_wp * epsmch ) then
          ! check this first to avoid possible overflow
          alpha = 1.0_wp
       else
          alpha = abs(dot_product(w%d,w%y_sharp))/ dSks
          alpha = min(1.0_wp,alpha)
       end if
       hf(:)  = alpha * hf(:)

       ! update S_k (again, as in N&W, Section 10.2)

       ! hf = hf + (1/yts) (y# - Sk d)^T y:
       alpha = 1.0_wp/yts
       call dGER(n,n,alpha,w%ysharpSks,1,w%y,1,hf,n)
       ! hf = hf + (1/yts) y^T (y# - Sk d):
       call dGER(n,n,alpha,w%y,1,w%ysharpSks,1,hf,n)
       ! hf = hf - ((y# - Sk d)^T d)/((yts)**2)) * y y^T
       alpha = -dot_product(w%ysharpSks,w%d)/(yts**2)
       call dGER(n,n,alpha,w%y,1,w%y,1,hf,n)
100    Continue
     end subroutine rank_one_update


     subroutine update_trust_region_radius(rho,options,inform,w)
       Implicit None
       real(wp), intent(inout) :: rho ! ratio of actual to predicted reduction
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_workspace ), intent(inout) :: w
       Integer :: nrec
       Character(Len=80) :: rec(1)
       Logical :: prnt3
       prnt3 = buildmsg(5,.false.,options)
       nrec = 0

       select case(options%tr_update_strategy)
       case(1) ! default, step-function
          if (rho < options%eta_success_but_reduce) then
             ! unsuccessful....reduce Delta
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             If (prnt3) Write(rec(1), Fmt=3010) w%Delta
             nrec=1
          else if (rho < options%eta_very_successful) then
             ! doing ok...retain status quo
             If (prnt3) Write(rec(1), Fmt=3020) w%Delta
             nrec=1
          else if (rho < options%eta_too_successful ) then
             ! more than very successful -- increase delta
             select case(options%type_of_method)
             case(1)
                w%Delta = min(options%maximum_radius, &
                     options%radius_increase * w%norm_S_d )
                ! if we have a trust region method, then we
                ! increase based on ||d||, not on Delta, as there's
                ! no point increasing the radius if we're within the
                ! trust region
             case(2)
                w%Delta = min(options%maximum_radius, &
                     options%radius_increase * w%Delta )
                ! increase based on Delta to ensure the
                ! regularized case works too
              Case Default
                ! This should never be reached under normal usage
                inform%status = NLLS_ERROR_UNSUPPORTED_TYPE_METHOD
                Go To 100
             end select
             If (prnt3) Write(rec(1), Fmt=3030) w%Delta
             nrec=1
          else if (rho < HUGE(wp)) then
             If (prnt3) Write(rec(1), Fmt=3040) w%Delta
             nrec=1
             ! too successful...accept step, but don't change w%Delta
          else
             ! just incase (NaNs, rho too big, and the like...)
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             rho = -1.0_wp ! set to be negative, so that the logic works....
             If (prnt3) Write(rec(1), Fmt=3050) w%Delta
             nrec=1
          end if
       case(2) ! Continuous method
          ! Based on that proposed by Hans Bruun Nielsen, TR IMM-REP-1999-05
          ! http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf
          if (rho >= options%eta_too_successful) then
             ! too successful....accept step, but don't change w%Delta
             If (prnt3) Write(rec(1), Fmt=3040) w%Delta
             nrec=1
          else if (rho > options%eta_successful) then
             w%Delta = w%Delta * min(options%radius_increase, &
                  max(options%radius_reduce, &
                  1 - ( (options%radius_increase - 1) * ((1 - 2*rho)**w%tr_p)) ))
             w%tr_nu = options%radius_reduce
             If (prnt3) Write(rec(1), Fmt=3060) w%Delta
             nrec=1
          else if ( rho <= options%eta_successful ) then
             w%Delta = w%Delta * w%tr_nu
             w%tr_nu =  w%tr_nu * 0.5_wp
             If (prnt3) Write(rec(1), Fmt=3010) w%Delta
             nrec=1
          else
             ! just incase (NaNs and the like...)
             w%Delta = max( options%radius_reduce, options%radius_reduce_max) * w%Delta
             rho = -1.0_wp ! set to be negative, so that the logic works....
             If (prnt3) Write(rec(1), Fmt=3050) w%Delta
             nrec=1
          end if
       case default
          inform%status = NLLS_ERROR_BAD_TR_STRATEGY
          Go To 100
       end select

100    Continue

       If (nrec>0.and.prnt3) Then
         Call Printmsg(3,.False.,options,nrec,rec)
       End If

3010   FORMAT('Unsuccessful step -- decreasing Delta to', ES12.4)
3020   FORMAT('Successful step -- Delta staying at', ES12.4)
3030   FORMAT('Very successful step -- increasing Delta to', ES12.4)
3040   FORMAT('Step too successful -- Delta staying at', ES12.4)
3050   FORMAT('NaN encountered -- reduced Delta to', ES12.4)
3060   FORMAT('Changing Delta to ', ES12.4)
     end subroutine update_trust_region_radius

     subroutine test_convergence(normF,normJF,normF0,normJF0,norm_2_d,options,inform)
       Implicit None
       real(wp), intent(in) :: normF, normJf, normF0, normJF0, norm_2_d
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       Character(Len=8), Parameter :: labels(3) = (/'||f||   ', 'gradient', 'step    '/)
       Integer :: nlabel
       Character(Len=80) :: rec(1)

       if ( normF <= max(options%stop_f_absolute, &
            options%stop_f_relative * normF0) ) then
          inform%convergence_normf = 1
          nlabel = 1
          Go To 100
       end if

       if ( (normJF/normF) <= max(options%stop_g_absolute, &
            options%stop_g_relative * (normJF0/normF0)) ) then
          inform%convergence_normg = 1
          nlabel = 2
          Go To 100
       end if
       if ( norm_2_d < options%stop_s ) then
          inform%convergence_norms = 1
          nlabel = 3
       end if

100   Continue
!     Pretty print results
      if ((inform%convergence_normf == 1 .Or. inform%convergence_normg == 1   &
         .Or. inform%convergence_norms == 1).And.buildmsg(5,.false.,options)) Then
        write(rec(1),Fmt=99999) 'Converged (',Trim(labels(nlabel)),   &
          ' test) at iteration', inform%iter
        Call Printmsg(5,.false.,options,1,rec)
99999 Format (3A,1X,I0)
      End If
     end subroutine test_convergence

     subroutine mult_J(J,n,m,x,Jx,options)
       Implicit None
       real(wp), intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       real(wp), intent(out) :: Jx(*)
       type(nlls_options), optional :: options

       real(wp) :: alpha, beta

       Jx(1:m) = 1.0_wp
       alpha = 1.0_wp
       beta  = 0.0_wp

       ! Avoid short-circuit evaluation
!      if ( present(options) .and. (.not. options%Fortran_Jacobian) ) then
       if (present(options)) Then
         If (.not. options%Fortran_Jacobian) then
           ! Jacobian held in row major format...
           call dgemv('T',n,m,alpha,J,n,x,1,beta,Jx,1)
         Else
           call dgemv('N',m,n,alpha,J,m,x,1,beta,Jx,1)
         End If
       else
          call dgemv('N',m,n,alpha,J,m,x,1,beta,Jx,1)
       end if
     end subroutine mult_J

     subroutine mult_Jt(J,n,m,x,Jtx,options)
       Implicit None
       double precision, intent(in) :: J(*), x(*)
       integer, intent(in) :: n,m
       double precision, intent(out) :: Jtx(*)
       type( nlls_options), optional :: options

       double precision :: alpha, beta

       Jtx(1:n) = 1.0_wp
       alpha = 1.0_wp
       beta  = 0.0_wp

       ! Avoid short-circuit evaluation
!      if (present(options).and. (.not. options%Fortran_Jacobian)) then
       if (present(options)) Then
         If (.not. options%Fortran_Jacobian) then
           ! Jacobian held in row major format...
           call dgemv('N',n,m,alpha,J,n,x,1,beta,Jtx,1)
         Else
           call dgemv('T',m,n,alpha,J,m,x,1,beta,Jtx,1)
         End If
       else
          call dgemv('T',m,n,alpha,J,m,x,1,beta,Jtx,1)
       end if
     end subroutine mult_Jt

     subroutine scale_J_by_weights(J,n,m,weights,options)
       Implicit None
       real(wp), intent(inout) :: J(*)
       real(wp), intent(in) :: weights(*)
       integer, intent(in) :: n,m
       type( nlls_options ) :: options

       integer :: i
       ! set J -> WJ
       if (options%Fortran_Jacobian) then
          do i = 1, n
             J( (i-1)*m + 1 : i*m) = weights(1:m)*J( (i-1)*m + 1 : i*m)
          end do
       else
          do i = 1, m
             J( (i-1)*n + 1 : i*n) = weights(i)*J( (i-1)*n + 1 : i*n)
          end do
       end if
     end subroutine scale_J_by_weights

!     subroutine add_matrices3(A,B,n,C)
!       Implicit None
!       real(wp), intent(in) :: A(*), B(*)
!       integer, intent(in) :: n
!       real(wp), intent(InOut) :: C(*)
!
!       C(1:n) = A(1:n) + B(1:n)
!     end subroutine add_matrices3

     subroutine add_matrices(A,B,n)
       Implicit None
       real(wp), intent(InOut) :: A(*)
       real(wp), intent(in) :: B(*)
       integer, intent(in) :: n
       A(1:n) = A(1:n) + B(1:n)
     end subroutine add_matrices

     subroutine get_element_of_matrix(J,m,ii,jj,Jij)
       Implicit None
       real(wp), intent(in) :: J(*)
       integer, intent(in) :: m
       integer, intent(in) :: ii,jj
       real(wp), intent(out) :: Jij

       ! return the (ii,jj)th entry of a matrix

       ! J held by columns....
       Jij = J(ii + (jj-1)*m)
     end subroutine get_element_of_matrix

     subroutine solve_spd(A,b,LtL,x,n,inform)
       Implicit None
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: LtL(:,:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform), intent(inout) :: inform
       inform%status = 0
       ! wrapper for the lapack subroutine dposv
       LtL(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dposv('L', n, 1, LtL, n, x, n, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dposv'
       end if
     end subroutine solve_spd

     subroutine solve_spd_nocopy(A,b,x,n,inform)
       Implicit None
       REAL(wp), intent(inout) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform), intent(inout) :: inform
       inform%status = 0
       ! wrapper for the lapack subroutine dposv
       x(1:n) = b(1:n)
       call dposv('L', n, 1, A, n, x, n, inform%external_return)
       if (inform%external_return .ne. 0) then
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dposv'
       end if
     end subroutine solve_spd_nocopy

     subroutine solve_general(A,b,x,n,inform,w)
       Implicit None
       REAL(wp), intent(in) :: A(:,:)
       REAL(wp), intent(in) :: b(:)
       REAL(wp), intent(out) :: x(:)
       integer, intent(in) :: n
       type( nlls_inform ), intent(inout) :: inform
       type( solve_general_work ) :: w
       ! wrapper for the lapack subroutine dposv
       ! NOTE: A would be destroyed

       If (.not. w%allocated) Then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       End If

       w%A(1:n,1:n) = A(1:n,1:n)
       x(1:n) = b(1:n)
       call dgesv( n, 1, w%A, n, w%ipiv, x, n, inform%external_return)
       if (inform%external_return .ne. 0 ) then
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dgesv'
       end if

100    continue
     end subroutine solve_general

     subroutine matrix_norm(x,A,norm_A_x)
       Implicit None
       REAL(wp), intent(in) :: A(:,:), x(:)
       REAL(wp), intent(out) :: norm_A_x
       ! Calculates norm_A_x = ||x||_A = sqrt(x'*A*x)
       norm_A_x = sqrt(dot_product(x,matmul(A,x)))
     end subroutine matrix_norm

     subroutine matmult_inner(J,n,m,A,options)
       Implicit None
       integer, intent(in) :: n,m
       real(wp), intent(in) :: J(*)
       real(wp), intent(out) :: A(n,n)
       type( nlls_options ), intent(in), optional :: options

       ! Takes an m x n matrix J and forms the
       ! n x n matrix A given by
       ! A = J' * J

       ! Avoid short-circuit evaluation
!      if ( present(options) .and. (.not. options%Fortran_Jacobian) ) then
       if ( present(options) ) Then
         If (.not. options%Fortran_Jacobian) then
           ! c format
           call dgemm('N','T',n, n, m, 1.0_wp, J, n, J, n, 0.0_wp, A, n)
         Else
           call dgemm('T','N',n, n, m, 1.0_wp, J, m, J, m, 0.0_wp, A, n)
         End If
       Else
          call dgemm('T','N',n, n, m, 1.0_wp, J, m, J, m, 0.0_wp, A, n)
       End If
     end subroutine matmult_inner
     subroutine matmult_outer(J,n,m,A)
       Implicit None
       integer, intent(in) :: n,m
       real(wp), intent(in) :: J(*)
       real(wp), intent(out) :: A(m,m)

       ! Takes an m x n matrix J and forms the
       ! m x m matrix A given by
       ! A = J * J'

       call dgemm('N','T',m, m, n, 1.0_wp, J, m, J, m, 0.0_wp, A, m)
     end subroutine matmult_outer

     subroutine outer_product(x,n,xxt)
       Implicit None
       real(wp), intent(in) :: x(:)
       integer, intent(in) :: n
       real(wp), intent(out) :: xxt(:,:)

       ! Takes an n vector x and forms the
       ! n x n matrix xtx given by
       ! xxt = x * x'

       xxt(1:n,1:n) = 0.0_wp
       call dger(n, n, 1.0_wp, x, 1, x, 1, xxt, n)
     end subroutine outer_product

     subroutine all_eig_symm(A,n,ew,ev,w,inform)
       Implicit None
       ! calculate all the eigenvalues of A (symmetric)
       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew(:), ev(:,:)
       type( all_eig_symm_work ) :: w
       type( nlls_inform ), intent(inout) :: inform
       integer :: lwork

       If (.not. w%allocated) Then
          inform%status = NLLS_ERROR_WORKSPACE_ERROR
          Go To 100
        End If

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
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dsyev'
       end if
100   continue
     end subroutine all_eig_symm

     subroutine min_eig_symm(A,n,ew,ev,options,inform,w)
       Implicit None
       ! calculate the leftmost eigenvalue of A
       real(wp), intent(in) :: A(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew, ev(:)
       type( nlls_inform ), intent(inout) :: inform
       type( nlls_options ), INTENT( IN ) :: options
       type( min_eig_symm_work ) :: w

       real(wp) :: tol, dlamch
       integer :: lwork, eigsout, minindex(1)

       If ( .not. w%allocated ) Then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       End If

       tol = 2.0_wp*dlamch('S')

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
             inform%status = NLLS_ERROR_FROM_EXTERNAL
             inform%external_name = 'lapack_dsyev'
             Go To 100
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
               1.0_wp, 1.0_wp, & ! not used for RANGE = 'I'
               1, 1, & ! only find the first eigenpair
               tol, & ! abstol for the eigensolver
               eigsout, & ! total number of eigs found
               w%ew, ev, & ! the eigenvalue and eigenvector
               n, & ! ldz (the eigenvector array)
               w%work, lwork, w%iwork, &  ! workspace
               w%ifail, & ! array containing indicies of non-converging ews
               inform%external_return)
          if (inform%external_return .ne. 0) then
             inform%status = NLLS_ERROR_FROM_EXTERNAL
             inform%external_name = 'lapack_dsyevx'
             Go To 100
          end if
          ew = w%ew(1)
       end if

100   continue
     end subroutine min_eig_symm

     subroutine max_eig(A,B,n,ew,ev,nullevs,options,inform,w)
       Implicit None
       real(wp), intent(inout) :: A(:,:), B(:,:)
       integer, intent(in) :: n
       real(wp), intent(out) :: ew, ev(:)
       real(wp), intent(inout), allocatable :: nullevs(:,:)
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( max_eig_work ) :: w

       integer :: lwork, maxindex(1), no_null, halfn
       real(wp):: tau
       integer :: i, ierr_dummy

       if (.not. w%allocated) then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       end if

       ! Find the max eigenvalue/vector of the generalized eigenproblem
       !     A * y = lam * B * y
       ! further, if ||y(1:n/2)|| \approx 0, find and return the
       ! eigenvectors y(n/2+1:n) associated with this

       ! check that n is even (important for hard case -- see below)
       if (modulo(n,2).ne.0) then
         inform%status = NLLS_ERROR_AINT_EIG_ODD
         Go To 100
       End If

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
          inform%status = NLLS_ERROR_FROM_EXTERNAL
          inform%external_name = 'lapack_dggev'
          Go To 100
       end if

       ! now find the rightmost real eigenvalue
       w%vecisreal = .true.
       where ( abs(w%alphaI) > 1.0e-8_wp ) w%vecisreal = .false.
       if (.not. any(w%vecisreal)) then
         inform%status = NLLS_ERROR_AINT_EIG_IMAG ! Eigs imaginary error
         goto 100
       end if
       w%ew_array(:) = w%alphaR(:)/w%beta(:)
       maxindex = maxloc(w%ew_array,mask=w%vecisreal)

       tau = 1.0e-4_wp ! todo -- pass this through from above...
       ! note n/2 always even -- validated by test on entry
       if (norm2( w%vr(1:halfn,maxindex(1)) ) < tau) then
          ! hard case
          ! let's find which ev's are null...
          w%nullindex = 0
          no_null = 0
          do i = 1,n
             if (norm2( w%vr(1:halfn,i)) < 1.0e-4_wp ) then
                no_null = no_null + 1
                w%nullindex(no_null) = i
             end if
          end do

!         allocate(nullevs(halfn,no_null))
          If (allocated(nullevs)) Then
            ! increase the size of the allocated array only if we need to
            if (no_null > size(nullevs,2)) then
              deallocate( nullevs, stat=ierr_dummy )
              allocate( nullevs(halfn,no_null) , stat = inform%alloc_status)
              if (inform%alloc_status /= 0) then
                inform%status = NLLS_ERROR_ALLOCATION
                inform%bad_alloc = "max_eig"
                goto 100
              End if
            end if
          Else ! size(nullevs,2) is considered 0
            ! Allocate space
            allocate( nullevs(halfn,no_null) , stat = inform%alloc_status)
            if (inform%alloc_status /= 0) then
              inform%status = NLLS_ERROR_ALLOCATION
              inform%bad_alloc = "max_eig"
              goto 100
            end if
          End If
          nullevs(1:halfn,1:no_null) = w%vr(halfn+1 : n,w%nullindex(1:no_null))
       end if

       ew = w%alphaR(maxindex(1))/w%beta(maxindex(1))
       ev(:) = w%vr(:,maxindex(1))

100   continue
     end subroutine max_eig

     subroutine shift_matrix(A,sigma,AplusSigma,n)
       Implicit None
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


     ! routines needed for the Newton tensor model
     subroutine solve_newton_tensor(J, f, eval_HF, X, n, m, Delta, &
                                    num_successful_steps, &
                                    d, md, params, options, inform, &
                                    w, tenJ, inner_workspace)
       Implicit None
       integer, intent(in)   :: n,m
       real(wp), intent(in)  :: f(:), J(:)
       real(wp) , intent(in) :: X(:), Delta
       integer, intent(in) :: num_successful_steps
       real(wp), intent(out) :: d(:)
       real(wp), intent(out) :: md
       procedure( eval_hf_type ) :: eval_HF
       class( params_base_type ), target :: params
       type( nlls_options ), intent(in) :: options
       type( nlls_inform ), intent(inout) :: inform
       type( solve_newton_tensor_work ) :: w
       type( tenJ_type ), intent(InOut), target :: tenJ
       type( NLLS_workspace ), Intent(InOut) :: inner_workspace

       type( nlls_inform ) :: tensor_inform
       integer :: i

       ! We need to solve the problem
       !   min 1/2 \sum_{i=1}^m t_{ik}^2(s) + 1/p \sigma_k ||s||^p_p
       ! where
       !   t_{ik}(s) := r_i(x_k) + s' g_i(x_k) + 1/2 s' B_ik s
       ! and where B_ik is a symmetric approx to Hi(x), the Hessian of r_i(x).
       ! To do this, we call ral_nlls recursively.
       ! First, we need to set up the eval_r/J/Hf functions needed here.

       If (.not. w%allocated) Then
         inform%status = NLLS_ERROR_WORKSPACE_ERROR
         goto 100
       End If

       ! save to params
       w%tparams%f(1:m) = f(1:m)
       w%tparams%Delta = Delta
       w%tparams%J(1:n*m) = J(1:n*m)
       w%tparams%X(1:n) = X(1:n)
       w%tparams%eval_HF => eval_HF
       w%tparams%parent_params => params
       w%tparams%tenJ => tenJ

       if (.not. w%tparams%eval_hp_provided) then
          ! let's get all the Hi's...
          do i = 1,m
             call get_Hi(n, m, X, params, i, w%tparams%Hi(:,:,i), eval_HF, inform)
             If (inform%status/=0) Then
               Go To 100
             End If
          end do
       end if

       d(1:n) = 0.0_wp

       ! send to ral_nlls_iterate to solve the subproblem recursively
       select case (options%inner_method)
       case (1) ! send in a base regularization parameter
          w%tensor_options%base_regularization = 1.0_wp / Delta
       case (2)  ! this uses p = 3, and solves a nlls problem in more unknowns explicitly
          w%tparams%p = 3.0_wp
          w%tparams%extra = 2
          d(1:n) = 1.0e-12_wp ! Hessian not defined at 0 if p /= 2, so set 'small'
       case (3)  ! this uses p = 2, and solves implicitly
          w%tensor_options%regularization_term = 1.0_wp / Delta
!          write(*,*) 'regularization_term = ', w%tensor_options%regularization_term
          w%m_in = m
       case (4) ! this uses p = 3, and solves implicitly
          d(1:n) = 1.0e-12_wp ! Hessian not defined at 0 if p /= 2, so set 'small'
          w%tensor_options%regularization_term = 1.0_wp / Delta
          w%m_in = m
       case (5) ! this uses p = 2, and solves a nlls problem in more unknowns explicitly
          w%tparams%p = 2.0_wp
          w%tparams%extra = 1
       end select

       tensor_inform%inner_iter = inform%iter + inform%inner_iter

       do i = 1, w%tensor_options%maxit
          call nlls_iterate(n,w%m_in,d, &
               inner_workspace, &
               evaltensor_f, evaltensor_J, evaltensor_HF, &
               w%tparams, &
               tensor_inform, w%tensor_options )
          if (tensor_inform%status /= 0) then
             ! there's an error : exit
             exit
          elseif ( (tensor_inform%convergence_normf == 1) &
               .or.(tensor_inform%convergence_normg == 1) &
               .or.(tensor_inform%convergence_norms == 1)) then
             ! we've converged!
             inform%inner_iter_success = .True.
             exit
          end if
       end do
       call nlls_finalize(inner_workspace,w%tensor_options)
       inform%inner_iter = inform%inner_iter + tensor_inform%iter

        ! now we need to evaluate the model at the new point
        w%tparams%extra = 0
        ! Does not fail
        call evaltensor_f(inform%external_return, n, m, d, &
             w%model_tensor, w%tparams)
        md = 0.5_wp * norm2( w%model_tensor(1:m) )**2
        ! + 0.5 * (1.0/Delta) * (norm2(d(1:n))**2)

100   continue
     end subroutine solve_newton_tensor

     subroutine evaltensor_f(status, n, m, s, f, params)
       Implicit None
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(out)   :: f
       class( params_base_type ), intent(inout) :: params

       ! Add default return value for status
       status = 0

       ! note:: tenJ is a global (but private to this module) derived type
       !        with components Hs and Js, which are real arrays
       ! The function we need to minimize is
       !  \sum_{i=1}^m t_ik(s) = 1/2 \sum_{i=1}^m (r_i(x_l) + s' g_i(x_k) + 1/2 s' B_ik s)^2

       select type(params)
       type is(tensor_params_type)
          ! Note: params%tenJ (tensor_params_type) points to NLLS_workspace%tenJ
          ! note: params%m contains 'm' from the original problem
          ! if we're passing in the reg. factor via the function/Jacobian, then
          ! 'm' here is m+n from the original problem

          ! First, calculate Js
          call mult_J(params%J(1:n*params%m),n,params%m,s,params%tenJ%Js)
          ! TODO: add options to allow calling from C

          ! Now, calculate s'Hs
          call calculate_sHs(n,m, s, params)

          ! put them all together for the first 1:m terms of f
          f(1:params%m) = params%f(1:params%m) + params%tenJ%Js(1:params%m) +  &
            0.5_wp*params%tenJ%stHs(1:params%m)

          if (params%extra == 1) then
             ! we're passing in the regularization via the function/Jacobian
             f(params%m + 1: params%m + n) = (1.0_wp/sqrt(params%Delta)) * s(1:n)
          elseif (params%extra == 2) then
             f(params%m + 1) = sqrt(2.0_wp/(params%Delta * params%p)) * &
                               (norm2(s(1:n))**(params%p/2.0_wp))
          end if
       end select
     end subroutine evaltensor_f

     subroutine calculate_sHs( n, m, s, params)
       Implicit None
       integer, intent(in) :: n, m
       real(wp), dimension(*), intent(in) :: s
       class( params_base_type ), intent(inout) :: params
       integer :: ii, status
       select type(params)
       type is(tensor_params_type)
          if (params%eval_hp_provided) then
             call params%eval_HP(status,n,params%m,params%x,s(1:n),params%tenJ%Hs,params%parent_params)
          else
             do ii = 1,params%m
                params%tenJ%H(1:n,1:n) = params%Hi(1:n,1:n,ii)
                params%tenJ%Hs(1:n,ii) = 0.0_wp
                call dgemv('N',n,n,1.0_wp,params%tenJ%H(1,1),n,s,1,0.0_wp,params%tenJ%Hs(1,ii),1)
             end do
          end if
          call dgemv('T',n,params%m,1.0_wp,params%tenJ%Hs(1,1),n,s,1,0.0_wp, params%tenJ%stHs(1),1)
       end select
     end subroutine calculate_sHs

     subroutine evaltensor_J(status, n, m, s, J, params)
       Implicit None
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(out)   :: J
       class( params_base_type ), intent(inout) :: params
       integer :: ii, jj
       ! Add default return value for status
       status = 0
       ! The function we need to return is
       !  g_i + H_i s
       ! where h_i and H_i are the gradient and hessians of the original problem

       select type(params)
       type is(tensor_params_type)
          ! note: params%m contains 'm' from the original problem
          ! if we're passing in the reg. factor via the function/Jacobian, then
          ! 'm' here is m+n from the original problem
          J(1:n*m) = 0.0_wp
          do jj = 1,n ! columns
             ! tenJ%Hs has been set by evaltensor_f, which is called first
             J( (jj-1)*m + 1 : (jj-1)*m + params%m) &
                  = params%J((jj-1)*params%m + 1 : jj*params%m) + params%tenJ%Hs(jj,1:params%m)
          end do
          if (params%extra == 1) then
             ! we're passing in the regularization via the function/Jacobian
             do ii = 1,n ! loop over the columns...
                J(m*(ii-1) + params%m + ii) = sqrt(1.0_wp/params%Delta)
             end do
          elseif (params%extra == 2) then
             do ii = 1, n ! loop over the columns....
                J(m*ii) = sqrt( (params%p)/(2.0_wp * params%Delta) ) * &
                                      (norm2(s(1:n))**( (params%p/2.0_wp) - 2.0_wp)) * &
                                      s(ii)
             end do
          end if
       end select
     end subroutine evaltensor_J

     subroutine evaltensor_HF(status, n, m, s, f, HF, params)
       Implicit None
       integer, intent(out) :: status
       integer, intent(in)  :: n
       integer, intent(in)  :: m
       real(wp), dimension(*), intent(in)    :: s
       real(wp), dimension(*), intent(in)   :: f
       real(wp), dimension(*), intent(out) :: HF
       class( params_base_type ), intent(inout) :: params
!!$    integer :: ii
       integer :: jj, kk
       real(wp) :: normx, hf_local
       status = 0
       select type(params)
       type is (tensor_params_type)
          params%tenJ%H(1:n,1:n) = 0.0_wp
          call params%eval_HF(status,n,params%m,params%x,f(1:m),HF(1:n**2),params%parent_params)
!!$          do ii = 1,params%m
!!$             params%tenJ%H(1:n,1:n) = params%tenJ%H(1:n,1:n) + f(ii)*params%Hi(1:n,1:n,ii)
!!$          end do
!!$
!!$          HF(1:n**2) = reshape(params%tenJ%H(1:n,1:n), (/n**2/))
          if (params%extra == 2) then
             normx = norm2(s(1:n))
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
       Implicit None
       integer, intent(in) :: n, m
       real(wp), intent(in) :: X(:)
       class( params_base_type ) :: params
       integer, intent(in) :: i
       real(wp), intent(out) :: Hi(:,:)
       procedure( eval_hf_type ) :: eval_HF
       type( nlls_inform ), intent( inout ) :: inform
       real( wp ), dimension( m ), intent( in ), optional :: weights
       real( wp ) :: ei( m )

       ei = 0.0_wp
       ei(i) = 1.0_wp

       if ( present(weights) ) then
          call eval_HF(inform%external_return, n, m, X, weights(1:m)*ei, Hi, params)
       else
          call eval_HF(inform%external_return, n, m, X, ei, Hi, params)
       end if
       If (inform%external_return/=0) Then
         inform%status = NLLS_ERROR_FROM_EXTERNAL
         inform%external_name = "eval_HF"
       End iF
     end subroutine get_Hi

     Subroutine check_options(opt, inform)
       ! Check some option values before enteting the solver
       Implicit None
       Type(NLLS_options), Intent(In) :: opt
       Type(NLLS_inform), Intent(InOut) :: inform
       Continue
       ! Routine to check for option value combinations
       If (opt%type_of_method < 1 .Or. Opt%type_of_method > 2) Then
         inform%status = NLLS_ERROR_UNSUPPORTED_TYPE_METHOD
       ElseIf (opt%print_level < 0 .Or. Opt%print_level > 5 ) Then
         inform%status = NLLS_ERROR_PRINT_LEVEL
!      ElseIf (opt%model==4.And.opt%nlls_method==3.And.opt%type_of_method==2) Then
!        If (opt%reg_order/=2.0_wp.And.opt%reg_order>0) Then
!         This specific option combination is not yet implemented!
!          inform%status = NLLS_ERROR_NOT_IMPLEMENTED
!        End If
!      ElseIf
!        ...
       End If

     End Subroutine

   end module ral_nlls_internal

