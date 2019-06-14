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

        nrec = 1
        Write (rec(nrec),Fmt=99995)
        nrec = nrec + 1
!       Begin OPTLIST.AWK contents
        Write(adj,Fmt=99999) "out"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%out
        nrec = nrec + 1
        Write(adj,Fmt=99999) "print_level"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%print_level
        nrec = nrec + 1
        Write(adj,Fmt=99999) "print_options"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%print_options
        nrec = nrec + 1
        Write(adj,Fmt=99999) "print_header"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%print_header
        nrec = nrec + 1
        Write(adj,Fmt=99999) "maxit"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%maxit
        nrec = nrec + 1
        Write(adj,Fmt=99999) "model"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%model
        nrec = nrec + 1
        Write(adj,Fmt=99999) "type_of_method"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%type_of_method
        nrec = nrec + 1
        Write(adj,Fmt=99999) "nlls_method"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%nlls_method
        nrec = nrec + 1
        Write(adj,Fmt=99999) "lls_solver"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%lls_solver
        nrec = nrec + 1
        Write(adj,Fmt=99999) "stop_g_absolute"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%stop_g_absolute
        nrec = nrec + 1
        Write(adj,Fmt=99999) "stop_g_relative"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%stop_g_relative
        nrec = nrec + 1
        Write(adj,Fmt=99999) "stop_f_absolute"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%stop_f_absolute
        nrec = nrec + 1
        Write(adj,Fmt=99999) "stop_f_relative"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%stop_f_relative
        nrec = nrec + 1
        Write(adj,Fmt=99999) "stop_s"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%stop_s
        nrec = nrec + 1
        Write(adj,Fmt=99999) "relative_tr_radius"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%relative_tr_radius
        nrec = nrec + 1
        Write(adj,Fmt=99999) "initial_radius_scale"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%initial_radius_scale
        nrec = nrec + 1
        Write(adj,Fmt=99999) "initial_radius"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%initial_radius
        nrec = nrec + 1
        Write(adj,Fmt=99999) "base_regularization"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%base_regularization
        nrec = nrec + 1
        Write(adj,Fmt=99999) "regularization"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%regularization
        nrec = nrec + 1
        Write(adj,Fmt=99999) "regularization_term"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%regularization_term
        nrec = nrec + 1
        Write(adj,Fmt=99999) "regularization_power"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%regularization_power
        nrec = nrec + 1
        Write(adj,Fmt=99999) "maximum_radius"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%maximum_radius
        nrec = nrec + 1
        Write(adj,Fmt=99999) "eta_successful"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%eta_successful
        nrec = nrec + 1
        Write(adj,Fmt=99999) "eta_success_but_reduce"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%eta_success_but_reduce
        nrec = nrec + 1
        Write(adj,Fmt=99999) "eta_very_successful"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%eta_very_successful
        nrec = nrec + 1
        Write(adj,Fmt=99999) "eta_too_successful"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%eta_too_successful
        nrec = nrec + 1
        Write(adj,Fmt=99999) "radius_increase"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%radius_increase
        nrec = nrec + 1
        Write(adj,Fmt=99999) "radius_reduce"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%radius_reduce
        nrec = nrec + 1
        Write(adj,Fmt=99999) "radius_reduce_max"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%radius_reduce_max
        nrec = nrec + 1
        Write(adj,Fmt=99999) "tr_update_strategy"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%tr_update_strategy
        nrec = nrec + 1
        Write(adj,Fmt=99999) "hybrid_switch"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%hybrid_switch
        nrec = nrec + 1
        Write(adj,Fmt=99999) "exact_second_derivatives"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%exact_second_derivatives
        nrec = nrec + 1
        Write(adj,Fmt=99999) "subproblem_eig_fact"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%subproblem_eig_fact
        nrec = nrec + 1
        Write(adj,Fmt=99999) "use_ews_subproblem"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%use_ews_subproblem
        nrec = nrec + 1
        Write(adj,Fmt=99999) "force_min_eig_symm"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%force_min_eig_symm
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%scale
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale_max"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%scale_max
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale_min"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%scale_min
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale_trim_min"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%scale_trim_min
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale_trim_max"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%scale_trim_max
        nrec = nrec + 1
        Write(adj,Fmt=99999) "scale_require_increase"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%scale_require_increase
        nrec = nrec + 1
        Write(adj,Fmt=99999) "setup_workspaces"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%setup_workspaces
        nrec = nrec + 1
        Write(adj,Fmt=99999) "remove_workspaces"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%remove_workspaces
        nrec = nrec + 1
        Write(adj,Fmt=99999) "more_sorensen_maxits"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%more_sorensen_maxits
        nrec = nrec + 1
        Write(adj,Fmt=99999) "more_sorensen_shift"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%more_sorensen_shift
        nrec = nrec + 1
        Write(adj,Fmt=99999) "more_sorensen_tiny"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%more_sorensen_tiny
        nrec = nrec + 1
        Write(adj,Fmt=99999) "more_sorensen_tol"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%more_sorensen_tol
        nrec = nrec + 1 
        Write(adj,Fmt=99999) "hybrid_tol"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%hybrid_tol
        nrec = nrec + 1
        Write(adj,Fmt=99999) "hybrid_switch_its"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%hybrid_switch_its
        nrec = nrec + 1
        Write(adj,Fmt=99999) "reg_order"
        Write(rec(nrec),Fmt=99998) Adjustl(adj), options%reg_order
        nrec = nrec + 1
        Write(adj,Fmt=99999) "inner_method"
        Write(rec(nrec),Fmt=99997) Adjustl(adj), options%inner_method
        nrec = nrec + 1
        Write(adj,Fmt=99999) "output_progress_vectors"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%output_progress_vectors
        nrec = nrec + 1
        Write(adj,Fmt=99999) "update_lower_order"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%update_lower_order
        nrec = nrec + 1
        Write(adj,Fmt=99999) "fortran_jacobian"
        Write(rec(nrec),Fmt=99996) Adjustl(adj), options%fortran_jacobian
        nrec = nrec + 1
        Write (rec(nrec),Fmt=99994)

!       End OPTLIST.AWK contents

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
