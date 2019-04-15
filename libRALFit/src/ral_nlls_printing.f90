    Module ral_nlls_printing
      Use ral_nlls_workspaces, Only: wp, nlls_options, nlls_inform
      Implicit None
      Private
      Public                           :: buildmsg, printmsg, print_options,   &
                                          print_bye

    Contains

      Function buildmsg(lev,onlythis,options)
!       @RALNAG buildm==e04svzn
!       options is a replacement for printing(:)
!       options%out == printing(1)
!       options%print_level == printing(2)

!       .. Implicit None Statement ..
        Implicit None
!       .. Function Return Value ..
        Logical                        :: buildmsg
!       .. Scalar Arguments ..
        Integer, Intent (In)           :: lev
        Logical, Intent (In)           :: onlythis
!       .. Array Arguments ..
        Type (nlls_options), Intent (In) :: options
!       .. Executable Statements ..
        Continue

        If (onlythis) Then
          buildmsg = options%print_level == lev
        Else
          buildmsg = options%print_level >= lev
        End If

        Return
      End Function buildmsg

!     Example on usage of printing prototype
      Subroutine printmsg(lev,onlythis,options,nrec,rec)
!       @RALNAG printm==e04svyn
!       options is a replacement for printing(:)
!       Intrinsic prints for RALFit can be disabled with enable=.False.

        Implicit None
        Integer, Intent (In)           :: lev, nrec
        Logical, Intent (In)           :: onlythis
        Type (nlls_options), Intent (In) :: options
        Character (*), Intent (In)     :: rec(nrec)

        Integer                        :: i

        Continue

        If (options%out>=0) Then
          If (options%print_level==lev .Or. .Not. onlythis .And. options%      &
            print_level>lev) Then
            Do i = 1, nrec
              Write (options%out,Fmt=99999) trim(rec(i))
            End Do
          End If
        End If

99999   Format (A)
      End Subroutine printmsg

      Subroutine print_options(options)
        Implicit None
        Type (nlls_options), Intent (In) :: options

        Integer                        :: nrec
        Character (80)                 :: rec(80)
        Character (30)                 :: adj

        Write (rec(1),Fmt=99995)
!       Begin OPTLIST.AWK contents
        Write(adj,Fmt=99999) "out"
        Write(rec(2),Fmt=99997) Adjustl(adj), options%out
        Write(adj,Fmt=99999) "print_level"
        Write(rec(3),Fmt=99997) Adjustl(adj), options%print_level
        Write(adj,Fmt=99999) "print_options"
        Write(rec(4),Fmt=99996) Adjustl(adj), options%print_options
        Write(adj,Fmt=99999) "print_header"
        Write(rec(5),Fmt=99997) Adjustl(adj), options%print_header
        Write(adj,Fmt=99999) "maxit"
        Write(rec(6),Fmt=99997) Adjustl(adj), options%maxit
        Write(adj,Fmt=99999) "model"
        Write(rec(7),Fmt=99997) Adjustl(adj), options%model
        Write(adj,Fmt=99999) "type_of_method"
        Write(rec(8),Fmt=99997) Adjustl(adj), options%type_of_method
        Write(adj,Fmt=99999) "nlls_method"
        Write(rec(9),Fmt=99997) Adjustl(adj), options%nlls_method
        Write(adj,Fmt=99999) "lls_solver"
        Write(rec(10),Fmt=99997) Adjustl(adj), options%lls_solver
        Write(adj,Fmt=99999) "stop_g_absolute"
        Write(rec(11),Fmt=99998) Adjustl(adj), options%stop_g_absolute
        Write(adj,Fmt=99999) "stop_g_relative"
        Write(rec(12),Fmt=99998) Adjustl(adj), options%stop_g_relative
        Write(adj,Fmt=99999) "stop_f_absolute"
        Write(rec(13),Fmt=99998) Adjustl(adj), options%stop_f_absolute
        Write(adj,Fmt=99999) "stop_f_relative"
        Write(rec(14),Fmt=99998) Adjustl(adj), options%stop_f_relative
        Write(adj,Fmt=99999) "stop_s"
        Write(rec(15),Fmt=99998) Adjustl(adj), options%stop_s
        Write(adj,Fmt=99999) "relative_tr_radius"
        Write(rec(16),Fmt=99997) Adjustl(adj), options%relative_tr_radius
        Write(adj,Fmt=99999) "initial_radius_scale"
        Write(rec(17),Fmt=99998) Adjustl(adj), options%initial_radius_scale
        Write(adj,Fmt=99999) "initial_radius"
        Write(rec(18),Fmt=99998) Adjustl(adj), options%initial_radius
        Write(adj,Fmt=99999) "base_regularization"
        Write(rec(19),Fmt=99998) Adjustl(adj), options%base_regularization
        Write(adj,Fmt=99999) "regularization"
        Write(rec(20),Fmt=99997) Adjustl(adj), options%regularization
        Write(adj,Fmt=99999) "regularization_term"
        Write(rec(21),Fmt=99998) Adjustl(adj), options%regularization_term
        Write(adj,Fmt=99999) "regularization_power"
        Write(rec(22),Fmt=99998) Adjustl(adj), options%regularization_power
        Write(adj,Fmt=99999) "maximum_radius"
        Write(rec(23),Fmt=99998) Adjustl(adj), options%maximum_radius
        Write(adj,Fmt=99999) "eta_successful"
        Write(rec(24),Fmt=99998) Adjustl(adj), options%eta_successful
        Write(adj,Fmt=99999) "eta_success_but_reduce"
        Write(rec(25),Fmt=99998) Adjustl(adj), options%eta_success_but_reduce
        Write(adj,Fmt=99999) "eta_very_successful"
        Write(rec(26),Fmt=99998) Adjustl(adj), options%eta_very_successful
        Write(adj,Fmt=99999) "eta_too_successful"
        Write(rec(27),Fmt=99998) Adjustl(adj), options%eta_too_successful
        Write(adj,Fmt=99999) "radius_increase"
        Write(rec(28),Fmt=99998) Adjustl(adj), options%radius_increase
        Write(adj,Fmt=99999) "radius_reduce"
        Write(rec(29),Fmt=99998) Adjustl(adj), options%radius_reduce
        Write(adj,Fmt=99999) "radius_reduce_max"
        Write(rec(30),Fmt=99998) Adjustl(adj), options%radius_reduce_max
        Write(adj,Fmt=99999) "tr_update_strategy"
        Write(rec(31),Fmt=99997) Adjustl(adj), options%tr_update_strategy
        Write(adj,Fmt=99999) "hybrid_switch"
        Write(rec(32),Fmt=99998) Adjustl(adj), options%hybrid_switch
        Write(adj,Fmt=99999) "exact_second_derivatives"
        Write(rec(33),Fmt=99996) Adjustl(adj), options%exact_second_derivatives
        Write(adj,Fmt=99999) "subproblem_eig_fact"
        Write(rec(34),Fmt=99996) Adjustl(adj), options%subproblem_eig_fact
        Write(adj,Fmt=99999) "use_ews_subproblem"
        Write(rec(35),Fmt=99996) Adjustl(adj), options%use_ews_subproblem
        Write(adj,Fmt=99999) "force_min_eig_symm"
        Write(rec(36),Fmt=99996) Adjustl(adj), options%force_min_eig_symm
        Write(adj,Fmt=99999) "scale"
        Write(rec(37),Fmt=99997) Adjustl(adj), options%scale
        Write(adj,Fmt=99999) "scale_max"
        Write(rec(38),Fmt=99998) Adjustl(adj), options%scale_max
        Write(adj,Fmt=99999) "scale_min"
        Write(rec(39),Fmt=99998) Adjustl(adj), options%scale_min
        Write(adj,Fmt=99999) "scale_trim_min"
        Write(rec(40),Fmt=99996) Adjustl(adj), options%scale_trim_min
        Write(adj,Fmt=99999) "scale_trim_max"
        Write(rec(41),Fmt=99996) Adjustl(adj), options%scale_trim_max
        Write(adj,Fmt=99999) "scale_require_increase"
        Write(rec(42),Fmt=99996) Adjustl(adj), options%scale_require_increase
        Write(adj,Fmt=99999) "calculate_svd_j"
        Write(rec(43),Fmt=99996) Adjustl(adj), options%calculate_svd_j
        Write(adj,Fmt=99999) "setup_workspaces"
        Write(rec(44),Fmt=99996) Adjustl(adj), options%setup_workspaces
        Write(adj,Fmt=99999) "remove_workspaces"
        Write(rec(45),Fmt=99996) Adjustl(adj), options%remove_workspaces
        Write(adj,Fmt=99999) "more_sorensen_maxits"
        Write(rec(46),Fmt=99997) Adjustl(adj), options%more_sorensen_maxits
        Write(adj,Fmt=99999) "more_sorensen_shift"
        Write(rec(47),Fmt=99998) Adjustl(adj), options%more_sorensen_shift
        Write(adj,Fmt=99999) "more_sorensen_tiny"
        Write(rec(48),Fmt=99998) Adjustl(adj), options%more_sorensen_tiny
        Write(adj,Fmt=99999) "more_sorensen_tol"
        Write(rec(49),Fmt=99998) Adjustl(adj), options%more_sorensen_tol
        Write(adj,Fmt=99999) "hybrid_tol"
        Write(rec(50),Fmt=99998) Adjustl(adj), options%hybrid_tol
        Write(adj,Fmt=99999) "hybrid_switch_its"
        Write(rec(51),Fmt=99997) Adjustl(adj), options%hybrid_switch_its
        Write(adj,Fmt=99999) "reg_order"
        Write(rec(52),Fmt=99998) Adjustl(adj), options%reg_order
        Write(adj,Fmt=99999) "inner_method"
        Write(rec(53),Fmt=99997) Adjustl(adj), options%inner_method
        Write(adj,Fmt=99999) "output_progress_vectors"
        Write(rec(54),Fmt=99996) Adjustl(adj), options%output_progress_vectors
        Write(adj,Fmt=99999) "update_lower_order"
        Write(rec(55),Fmt=99996) Adjustl(adj), options%update_lower_order
        Write(adj,Fmt=99999) "fortran_jacobian"
        Write(rec(56),Fmt=99996) Adjustl(adj), options%fortran_jacobian
        Write(adj,Fmt=99999) "box_nfref_max"
        Write(rec(57),Fmt=99997) Adjustl(adj), options%box_nfref_max
        Write(adj,Fmt=99999) "box_ntrfail"
        Write(rec(58),Fmt=99997) Adjustl(adj), options%box_ntrfail
        Write(adj,Fmt=99999) "box_gamma"
        Write(rec(59),Fmt=99998) Adjustl(adj), options%box_gamma
        Write(adj,Fmt=99999) "box_decmin"
        Write(rec(60),Fmt=99998) Adjustl(adj), options%box_decmin
        Write(adj,Fmt=99999) "box_bigbnd"
        Write(rec(61),Fmt=99998) Adjustl(adj), options%box_bigbnd
        Write(adj,Fmt=99999) "box_wolfe_descent"
        Write(rec(62),Fmt=99998) Adjustl(adj), options%box_wolfe_descent
        Write(adj,Fmt=99999) "box_wolfe_curvature"
        Write(rec(63),Fmt=99998) Adjustl(adj), options%box_wolfe_curvature
        Write(adj,Fmt=99999) "box_kanzow_power"
        Write(rec(64),Fmt=99998) Adjustl(adj), options%box_kanzow_power
        Write(adj,Fmt=99999) "box_kanzow_descent"
        Write(rec(65),Fmt=99998) Adjustl(adj), options%box_kanzow_descent
        Write(adj,Fmt=99999) "box_quad_model_descent"
        Write(rec(66),Fmt=99998) Adjustl(adj), options%box_quad_model_descent
        Write(adj,Fmt=99999) "box_tr_test_step"
        Write(rec(67),Fmt=99996) Adjustl(adj), options%box_tr_test_step
        Write(adj,Fmt=99999) "box_wolfe_test_step"
        Write(rec(68),Fmt=99996) Adjustl(adj), options%box_wolfe_test_step
        Write(adj,Fmt=99999) "box_tau_max"
        Write(rec(69),Fmt=99998) Adjustl(adj), options%box_tau_max
        Write(adj,Fmt=99999) "box_max_ntrfail"
        Write(rec(70),Fmt=99997) Adjustl(adj), options%box_max_ntrfail
        Write(adj,Fmt=99999) "box_quad_match"
        Write(rec(71),Fmt=99997) Adjustl(adj), options%box_quad_match
        Write(adj,Fmt=99999) "box_alpha_scale"
        Write(rec(72),Fmt=99998) Adjustl(adj), options%box_alpha_scale
        Write(adj,Fmt=99999) "box_delta_scale"
        Write(rec(73),Fmt=99998) Adjustl(adj), options%box_delta_scale
        Write(adj,Fmt=99999) "box_tau_min"
        Write(rec(74),Fmt=99998) Adjustl(adj), options%box_tau_min
        Write(adj,Fmt=99999) "box_ls_step_maxit"
        Write(rec(75),Fmt=99997) Adjustl(adj), options%box_ls_step_maxit
        Write(adj,Fmt=99999) "box_linesearch_type"
        Write(rec(76),Fmt=99997) Adjustl(adj), options%box_linesearch_type
        Write (rec(77),Fmt=99994)
        nrec = 77
!       End OPTLIST.AWK contents

        Call printmsg(0,.False.,options,nrec,rec)

99999   Format (A30)
99998   Format (5X,A30,'=',8X,Es12.4e3)
99997   Format (5X,A30,'=',10X,I10)
99996   Format (5X,A30,'=',10X,L10)
99995   Format (1X,'Begin of Options')
99994   Format (1X,'End of Options')

      End Subroutine print_options


      Subroutine print_bye(options,inform,box)
!       Assumes print level > 0
        Implicit None
        Type (nlls_options), Intent (In) :: options
        Type (nlls_inform), Intent (In) :: inform
        Logical, Intent (In)           :: box

        Integer                        :: nrec, tot
        Character (90)                 :: rec(20)

        Continue

        If (buildmsg(2,.False.,options)) Then
          nrec = 1
          Write (rec(nrec),Fmt=99993)
          Call printmsg(2,.False.,options,nrec,rec)
        End If

        If (inform%status/=0) Then
          nrec = 1
          Write (rec(nrec),Fmt=99994) 'Error: ', inform%status
          nrec = nrec + 1
          Write (rec(nrec),Fmt=99993)
          nrec = nrec + 1
          Write (rec(nrec),Fmt=99992) trim(inform%error_message)
          nrec = nrec + 1
          Write (rec(nrec),Fmt=99991)
          Call printmsg(1,.False.,options,nrec,rec)
        Else
          If (buildmsg(1,.False.,options)) Then
            nrec = 1
            Write (rec(nrec),Fmt=99998)                                        &
              'converged, an optimal solution was found'
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99993)
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99997) 'Norm of error                 ',      &
              inform%obj
            Call printmsg(1,.False.,options,nrec,rec)
          End If
          If (buildmsg(2,.False.,options)) Then
            tot = max(1,inform%iter+inform%ls_step_iter+inform%pg_step_iter)
            nrec = 1
            If (box) Then
              Write (rec(nrec),Fmt=99997) 'Norm of projected gradient    ',    &
                inform%norm_g
              nrec = nrec + 1
              Write (rec(nrec),Fmt=99997) 'Norm of scaled proj. gradient ',    &
                inform%scaled_g
            Else
              Write (rec(nrec),Fmt=99997) 'Norm of gradient              ',    &
                inform%norm_g
              nrec = nrec + 1
              Write (rec(nrec),Fmt=99997) 'Norm of scaled gradient       ',    &
                inform%scaled_g
            End If
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99996) 'Iteration count               ', tot
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    Trust region step         ',      &
              inform%iter, 100.0_wp*inform%iter/tot
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    LS step                   ',      &
              inform%ls_step_iter, 100.0_wp*inform%ls_step_iter/tot
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    PG step                   ',      &
              inform%pg_step_iter, 100.0_wp*inform%pg_step_iter/tot
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99996) 'Function evaluations          ',      &
              inform%f_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    Trust region step         ',      &
              inform%f_eval - inform%f_eval_ls - inform%f_eval_pg,             &
              100.0_wp*(inform%f_eval-inform%f_eval_ls-inform%f_eval_pg)/      &
              inform%f_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    LS step                   ',      &
              inform%f_eval_ls, 100.0_wp*real(inform%f_eval_ls)/inform%f_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    PG step                   ',      &
              inform%f_eval_pg, 100.0_wp*real(inform%f_eval_pg)/inform%f_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99996) 'Gradient evaluations          ',      &
              inform%g_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    Trust region step         ',      &
              inform%g_eval - inform%g_eval_ls - inform%g_eval_pg,             &
              100.0_wp*(inform%g_eval-inform%g_eval_ls-inform%g_eval_pg)/      &
              inform%g_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    LS step                   ',      &
              inform%g_eval_ls, 100.0_wp*real(inform%g_eval_ls)/inform%g_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99995) '    PG step                   ',      &
              inform%g_eval_pg, 100.0_wp*real(inform%g_eval_pg)/inform%g_eval
            nrec = nrec + 1
            Write (rec(nrec),Fmt=99996) 'Hessian evaluations           ',      &
              inform%h_eval
            Call printmsg(2,.False.,options,nrec,rec)
          End If
        End If

99999   Format (A)
99998   Format (2X,'Status:',1X,A)
99997   Format (1X,A30,4X,Es12.5e2)
99996   Format (1X,A30,4X,I12)
99995   Format (1X,A30,4X,I12,1X,'(',F5.1,'%)')
99994   Format (1X,A,4X,I0)
99993   Format (1X,57('-'))

99992   Format (1X,'**',1X,A)
99991   Format (1X,'** ABNORMAL EXIT from RALFit.')
      End Subroutine print_bye

    End Module ral_nlls_printing
