    Module ral_nlls_printing
      Use ral_nlls_workspaces, Only: nlls_options, nlls_inform
      Implicit None
      Private
      Public                           :: buildmsg, printmsg, print_options, print_bye

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

        Write (adj,Fmt=99999) 'error'
        Write (rec(2),Fmt=99997) adjustl(adj), options%error
        Write (adj,Fmt=99999) 'out'
        Write (rec(3),Fmt=99997) adjustl(adj), options%out
        Write (adj,Fmt=99999) 'print_level'
        Write (rec(4),Fmt=99997) adjustl(adj), options%print_level
        Write (adj,Fmt=99999) 'maxit'
        Write (rec(5),Fmt=99997) adjustl(adj), options%maxit
        Write (adj,Fmt=99999) 'model'
        Write (rec(6),Fmt=99997) adjustl(adj), options%model
        Write (adj,Fmt=99999) 'type_of_method'
        Write (rec(7),Fmt=99997) adjustl(adj), options%type_of_method
        Write (adj,Fmt=99999) 'nlls_method'
        Write (rec(8),Fmt=99997) adjustl(adj), options%nlls_method
        Write (adj,Fmt=99999) 'lls_solver'
        Write (rec(9),Fmt=99997) adjustl(adj), options%lls_solver
        Write (adj,Fmt=99999) 'stop_g_absolute'
        Write (rec(10),Fmt=99998) adjustl(adj), options%stop_g_absolute
        Write (adj,Fmt=99999) 'stop_g_relative'
        Write (rec(11),Fmt=99998) adjustl(adj), options%stop_g_relative
        Write (adj,Fmt=99999) 'stop_f_absolute'
        Write (rec(12),Fmt=99998) adjustl(adj), options%stop_f_absolute
        Write (adj,Fmt=99999) 'stop_f_relative'
        Write (rec(13),Fmt=99998) adjustl(adj), options%stop_f_relative
        Write (adj,Fmt=99999) 'stop_s'
        Write (rec(14),Fmt=99998) adjustl(adj), options%stop_s
        Write (adj,Fmt=99999) 'relative_tr_radius'
        Write (rec(15),Fmt=99997) adjustl(adj), options%relative_tr_radius
        Write (adj,Fmt=99999) 'initial_radius_scale'
        Write (rec(16),Fmt=99998) adjustl(adj), options%initial_radius_scale
        Write (adj,Fmt=99999) 'initial_radius'
        Write (rec(17),Fmt=99998) adjustl(adj), options%initial_radius
        Write (adj,Fmt=99999) 'base_regularization'
        Write (rec(18),Fmt=99998) adjustl(adj), options%base_regularization
        Write (adj,Fmt=99999) 'regularization_term'
        Write (rec(19),Fmt=99998) adjustl(adj), options%regularization_term
        Write (adj,Fmt=99999) 'regularization_power'
        Write (rec(20),Fmt=99998) adjustl(adj), options%regularization_power
        Write (adj,Fmt=99999) 'maximum_radius'
        Write (rec(21),Fmt=99998) adjustl(adj), options%maximum_radius
        Write (adj,Fmt=99999) 'eta_successful'
        Write (rec(22),Fmt=99998) adjustl(adj), options%eta_successful
        Write (adj,Fmt=99999) 'eta_success_but_reduce'
        Write (rec(23),Fmt=99998) adjustl(adj), options%eta_success_but_reduce
        Write (adj,Fmt=99999) 'eta_very_successful'
        Write (rec(24),Fmt=99998) adjustl(adj), options%eta_very_successful
        Write (adj,Fmt=99999) 'eta_too_successful'
        Write (rec(25),Fmt=99998) adjustl(adj), options%eta_too_successful
        Write (adj,Fmt=99999) 'radius_increase'
        Write (rec(26),Fmt=99998) adjustl(adj), options%radius_increase
        Write (adj,Fmt=99999) 'radius_reduce'
        Write (rec(27),Fmt=99998) adjustl(adj), options%radius_reduce
        Write (adj,Fmt=99999) 'radius_reduce_max'
        Write (rec(28),Fmt=99998) adjustl(adj), options%radius_reduce_max
        Write (adj,Fmt=99999) 'hybrid_switch'
        Write (rec(29),Fmt=99998) adjustl(adj), options%hybrid_switch
        Write (adj,Fmt=99999) 'exact_second_derivatives'
        Write (rec(30),Fmt=99996) adjustl(adj), options%                       &
          exact_second_derivatives
        Write (adj,Fmt=99999) 'subproblem_eig_fact'
        Write (rec(31),Fmt=99996) adjustl(adj), options%subproblem_eig_fact
        Write (adj,Fmt=99999) 'use_ews_subproblem'
        Write (rec(32),Fmt=99996) adjustl(adj), options%use_ews_subproblem
        Write (adj,Fmt=99999) 'scale'
        Write (rec(33),Fmt=99997) adjustl(adj), options%scale
        Write (adj,Fmt=99999) 'scale_max'
        Write (rec(34),Fmt=99998) adjustl(adj), options%scale_max
        Write (adj,Fmt=99999) 'scale_min'
        Write (rec(35),Fmt=99998) adjustl(adj), options%scale_min
        Write (adj,Fmt=99999) 'scale_trim_min'
        Write (rec(36),Fmt=99996) adjustl(adj), options%scale_trim_min
        Write (adj,Fmt=99999) 'scale_trim_max'
        Write (rec(37),Fmt=99996) adjustl(adj), options%scale_trim_max
        Write (adj,Fmt=99999) 'scale_require_increase'
        Write (rec(38),Fmt=99996) adjustl(adj), options%scale_require_increase
        Write (adj,Fmt=99999) 'calculate_svd_J'
        Write (rec(39),Fmt=99996) adjustl(adj), options%calculate_svd_j
        Write (adj,Fmt=99999) 'setup_workspaces'
        Write (rec(40),Fmt=99996) adjustl(adj), options%setup_workspaces
        Write (adj,Fmt=99999) 'remove_workspaces'
        Write (rec(41),Fmt=99996) adjustl(adj), options%remove_workspaces
        Write (adj,Fmt=99999) 'more_sorensen_maxits'
        Write (rec(42),Fmt=99997) adjustl(adj), options%more_sorensen_maxits
        Write (adj,Fmt=99999) 'more_sorensen_shift'
        Write (rec(43),Fmt=99998) adjustl(adj), options%more_sorensen_shift
        Write (adj,Fmt=99999) 'more_sorensen_tiny'
        Write (rec(44),Fmt=99998) adjustl(adj), options%more_sorensen_tiny
        Write (adj,Fmt=99999) 'more_sorensen_tol'
        Write (rec(45),Fmt=99998) adjustl(adj), options%more_sorensen_tol
        Write (adj,Fmt=99999) 'hybrid_tol'
        Write (rec(46),Fmt=99998) adjustl(adj), options%hybrid_tol
        Write (adj,Fmt=99999) 'hybrid_switch_its'
        Write (rec(47),Fmt=99997) adjustl(adj), options%hybrid_switch_its
        Write (adj,Fmt=99999) 'reg_order'
        Write (rec(48),Fmt=99998) adjustl(adj), options%reg_order
        Write (adj,Fmt=99999) 'inner_method'
        Write (rec(49),Fmt=99997) adjustl(adj), options%inner_method
        Write (adj,Fmt=99999) 'output_progress_vectors'
        Write (rec(50),Fmt=99996) adjustl(adj), options%                       &
          output_progress_vectors
        Write (adj,Fmt=99999) 'update_lower_order'
        Write (rec(51),Fmt=99996) adjustl(adj), options%update_lower_order
        Write (adj,Fmt=99999) 'Fortran_Jacobian'
        Write (rec(52),Fmt=99996) adjustl(adj), options%fortran_jacobian

        Write (rec(53),Fmt=99994)

        nrec = 53

        Call printmsg(0,.False.,options,nrec,rec)

99999   Format (A30)
99998   Format (5X,A30,'=',8X,Es12.4e3)
99997   Format (5X,A30,'=',10X,I10)
99996   Format (5X,A30,'=',10X,L10)
99995   Format (1X,'Begin of Options')
99994   Format (1X,'End of Options')

      End Subroutine print_options


      Subroutine print_bye(options,inform)
!       Assumes print level > 0
        Implicit None
        Type (nlls_options), Intent (In) :: options
        Type (nlls_inform), Intent (In) :: inform

        Integer                        :: nrec
        Character (Len=80)                 :: rec(10)

        Continue

     If (buildmsg(2,.False.,options)) Then
          nrec = 1
          Write(rec(nrec), Fmt=60000)
         Call printmsg(2, .False., options, nrec, rec)
     End If

     If (inform%status /= 0) then
         nrec = 1
         Write(rec(nrec), Fmt=99990) 'Error: ', inform%status
         nrec = nrec + 1
         Write(rec(nrec), Fmt=60000)
         nrec = nrec + 1
         Write(rec(nrec), Fmt=5000) trim(inform%error_message)
         nrec = nrec + 1
         Write(rec(nrec), Fmt=5001) inform%status
         Call printmsg(1, .False., options, nrec, rec)
     Else
        If (buildmsg(1,.False.,options)) Then
          nrec = 1
          Write(rec(nrec), Fmt=99998) 'converged, an optimal solution was found'
          nrec = nrec + 1
          Write(rec(nrec), Fmt=60000)
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99997) 'Norm of error                 ', inform%obj
          Call printmsg(1,.False.,options,nrec,rec)
        End If
        If (buildmsg(2,.False.,options)) Then
          nrec = 1
          Write(rec(nrec), Fmt=99997) 'Norm of gradient              ', inform%norm_g
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99997) 'Norm of scaled gradient       ', inform%scaled_g
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99997) 'Step size                     ', inform%step
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99996) 'Iteration count               ', inform%iter
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99996) 'Function evaluations          ', inform%f_eval
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99996) 'Gradient evaluations          ', inform%g_eval
          nrec = nrec + 1
          Write(rec(nrec), Fmt=99996) 'Hessian evaluations           ', inform%h_eval
          Call printmsg(2,.False.,options,nrec,rec)
        End If
     End If

99999 Format(A)
99998 Format(2X,'Status:',1X,A)
99997 Format(1X,A30,4X,Es12.5e2)
99996 Format(1X,A30,4X,I12)
99990 Format(1X,A,4X,I0)
60000 Format(1X,53('-'))

5000 Format(1X,'**',1X,A)
5001 Format(1X,'** ABNORMAL EXIT from RALFit routine nlls_solve: ERROR =',I5)
      End Subroutine print_bye

    End Module ral_nlls_printing
