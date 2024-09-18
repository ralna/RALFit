! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors
!    may be used to endorse or promote products derived from this software without
!    specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
! OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! ral_nlls_fd :: Finite differences
module ral_nlls_fd

   Use ral_nlls_workspaces, Only : params_base_type, params_internal_type,     &
      nlls_options, nlls_inform, nlls_workspace, wp, box_type

   Implicit None

   Public :: eval_f_wrap, eval_j_wrap, eval_hf_wrap, eval_hp_wrap, check_jacobian
   ! Public interface to estimate the jacobian
   Public :: jacobian_calc, jacobian_setup, jacobian_free, jacobian_handle
   ! Public convinience default dummies for call-backs
   Public :: ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy

   Private

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

   ! handle
   Type :: jacobian_handle
      Logical :: setup = .False.
      Integer :: n, m
      Procedure(eval_f_type), NoPass, Pointer :: eval_f => Null()
      Class(params_base_type), Pointer :: params => Null()
      Real(Kind=wp), Pointer :: lower => Null(), upper => Null()
      Logical :: f_storage
      ! private elements for internal params
      Type(params_internal_type) :: iparams
      Type(nlls_inform) :: inform
      Type(NLLS_options) :: options
      Type(box_type) :: box
   End Type jacobian_handle

Contains

   ! Fortran dummy eval_J routine to request Finite-Differences Jacobian
   subroutine ral_nlls_eval_j_dummy(status, n, m, x, J, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(out) :: J
      class(params_base_type), intent(inout) :: params

      status = -45544554 ! Magic number to request FD
   end subroutine ral_nlls_eval_j_dummy

   ! Fortran dummy eval_HF routine
   subroutine ral_nlls_eval_hf_dummy(status, n, m, x, f, h, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(in)  :: f
      double precision, dimension(*), intent(out) :: h
      class(params_base_type), intent(inout) :: params

      status = -1023
   end subroutine ral_nlls_eval_hf_dummy

   Subroutine jacobian_setup(status, handle, n, m, x, eval_f, params, lower, upper, f_storage)
      Use ral_nlls_workspaces, Only : params_internal_type, wp, NLLS_options, &
         nlls_inform, setup_bounds_type, setup_iparams_type, NLLS_ERROR_ALLOCATION
      Implicit None
      Integer, Intent(out) :: status
      Type(jacobian_handle), Intent(Inout) :: handle
      Integer, Intent(in) :: n, m
      Real(Kind=wp), Dimension(n), Intent(In) :: x
      Procedure( eval_f_type ) :: eval_f
      Class(params_base_type), Intent(inout), Optional :: params
      Real(Kind=wp), Dimension(n), Intent(in), Target, Optional :: lower, upper
      Logical, Optional, Intent(in) :: f_storage

      Type(params_base_type) :: no_params
      Integer :: ierr_dummy

      Continue

      If (n<1 .Or. m <1) Then
         status = 1
         GoTo 100
      End If

      handle%n = n
      handle%m = m

      Call setup_bounds_type(handle%box,n,lower_bounds=lower, upper_bounds=upper, &
         options=handle%options, inform=handle%inform)
      status = handle%inform%status
      If (status /= 0) Then
         GoTo 100
      End If
      ! Free unrequired vectors
      If (allocated(handle%box%pdir)) deallocate(handle%box%pdir, Stat=ierr_dummy)
      If (allocated(handle%box%normFref)) deallocate(handle%box%normFref, Stat=ierr_dummy)
      If (allocated(handle%box%sk)) deallocate(handle%box%sk, Stat=ierr_dummy)
      If (allocated(handle%box%g)) deallocate(handle%box%g, Stat=ierr_dummy)

      ! Note x is not passed to setup iparams, it is assigned along with f in jacobian_calc
      If (Present(params)) Then
         Call setup_iparams_type(m, handle%iparams, params, eval_f,   &
            ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, handle%inform, handle%options, handle%box)
      Else
         Call setup_iparams_type(m, handle%iparams, no_params, eval_f,   &
            ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, handle%inform, handle%options, handle%box)
      End If

      ! Allocation must come AFTER call to setup_iparams_type()
      Allocate(handle%iparams%f_pert(m) ,Stat=status)
      If (status /= 0) Then
         status = NLLS_ERROR_ALLOCATION
         GoTo 100
      End If

      If (Present(f_storage)) Then
         handle%iparams%options%Fortran_Jacobian = f_storage
      End If

      handle%setup = .True.
      status = 0
      Return
100   Continue
      Call jacobian_free(handle)
   End subroutine jacobian_setup

   Subroutine jacobian_free(handle)
      Use ral_nlls_workspaces, Only : free_iparams_type, remove_workspace_bounds
      Implicit None
      Type(jacobian_handle), Intent(Inout) :: handle

      Continue

      handle%setup = .False.
      Call free_iparams_type(handle%iparams)
      Call remove_workspace_bounds(handle%box)
      handle%eval_f => Null()
      handle%params => Null()
      handle%lower => Null()
      handle%upper => Null()
   End Subroutine jacobian_free

   ! Assumes x and f are allocated and than f = f(x)
   subroutine jacobian_calc(status, handle, x, f, J, fd_step)
      use ral_nlls_workspaces, only : wp, NLLS_ERROR_WORKSPACE_ERROR
      implicit none
      integer, intent(out) :: status
      Type(jacobian_handle), Intent(Inout) :: handle
      Real(Kind=wp), dimension(*), intent(inout) :: x
      Real(Kind=wp), dimension(*), intent(in) :: f
      Real(Kind=wp), dimension(*), intent(out) :: J
      Real(Kind=wp), optional :: fd_step

      Continue

      If (.Not. handle%setup) Then
         status =  NLLS_ERROR_WORKSPACE_ERROR
         return
      End If

      If (present(fd_step)) Then
         handle%iparams%options%fd_step = fd_step
      End If

      Call fd_jacobian(status, handle%n, handle%m, x, f, J, handle%iparams)

      return
   end subroutine jacobian_calc

   ! Add an internal wrapper for eval_J (used to provide finite-differences)
   subroutine eval_j_wrap(status, n, m, x, J, params)
      use ral_nlls_workspaces, only : params_base_type, params_internal_type,  &
         nlls_options, nlls_inform, nlls_workspace
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(out) :: J
      class(params_base_type), intent(inout) :: params

      Continue

      ! Make sure we have the correct object for params
      Select Type( params )
       type is ( params_internal_type )
         if ( params%fd_type == 'N' ) then
            ! call exact Jacobian
            Call params%user%eval_J(status, n, m, x(1:n), J(1:n*m), params%user%params)
            if (status == -45544554) then
               ! Magic number from eval_j_dummy to indicate to use FD
               ! This branch should only be executed once
               ! Allocate space for f_pert
               Allocate(params%f_pert(m), Stat=status)
               if (status /= 0) return
               params%fd_type = 'Y'
            else
               ! User provided a custom eval_J, return user status and J values
               return
            end if
         end if
         if ( params%fd_type /= 'N' ) then
            ! Finite-differences
            ! Assumes that eval_F was called previously and that f(:) is valid
            Call fd_jacobian(status, n, m, params%x, params%f, J, params)
            return
         end if
      End Select
      status = -2024 ! Inform Jacobian could not be evaluated
   end subroutine eval_j_wrap


   ! transparent wrappers for eval_f, eval_h*
   subroutine eval_f_wrap(status, n, m, x, f, params)
      use ral_nlls_workspaces, only : params_base_type, params_internal_type
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(out) :: f
      class(params_base_type), intent(inout) :: params

      ! Make sure we have the correct object for params
      Select Type( params )
       type is ( params_internal_type)
         Call params%user%eval_f(status, n, m, x(1:n), f(1:m), params%user%params)
         return
      End Select
      status = -2023 ! Inform call-back could not be evaluated
   end subroutine eval_f_wrap

   subroutine eval_hf_wrap(status, n, m, x, f, h, params)
      use ral_nlls_workspaces, only : params_base_type, params_internal_type
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(in)  :: f
      double precision, dimension(*), intent(out) :: h
      class(params_base_type), intent(inout) :: params

      ! Make sure we have the correct object for params
      Select Type( params )
       type is ( params_internal_type)
         Call params%user%eval_hf(status, n, m, x(1:n), f(1:m), h(1:n*n), params%user%params)
         return
      End Select
      status = -2022 ! Inform call-back could not be evaluated
   end subroutine eval_hf_wrap

   subroutine eval_hp_wrap(status, n, m, x, y, hp, params)
      use ral_nlls_workspaces, only : params_base_type, params_internal_type
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(in)  :: y
      double precision, dimension(*), intent(out) :: hp
      class(params_base_type), intent(inout) :: params

      ! Make sure we have the correct object for params
      Select Type( params )
       type is ( params_internal_type)
         Call params%user%eval_hp(status, n, m, x(1:n), y(1:n), hp(1:n*m), params%user%params)
         return
      End Select
      status = -2021 ! Inform call-back could not be evaluated
   end subroutine eval_hp_wrap

   ! JACOBIAN DERIVATIVE CHECKER
   ! Compares the user call-back provided Jacobian agains a FD approximation
   ! The check is done is the caller is from nlls_solve only
   Subroutine check_jacobian(n, m, J, iparams)
      Use ral_nlls_workspaces, only: NLLS_ERROR_UNEXPECTED, &
         NLLS_ERROR_ALLOCATION, NLLS_ERROR_FROM_EXTERNAL, &
         NLLS_ERROR_BAD_JACOBIAN
      Use ral_nlls_printing, only: printmsg
      Implicit None
      Integer, Intent(In) :: n, m
      Real(Kind=wp), Dimension(:), Intent(In) :: J
      Type(params_internal_type), Intent(Inout) :: iparams

      Real(Kind=wp), Allocatable, Dimension(:) :: J_fd
      Integer :: ierr, ivar, jcon, iivar, jjcon, idx, idx_tran, prlvl, tcnt, skip
      Real(Kind=wp) :: perturbation, this_perturbation, relerr, relerr_tran, test_tol
      Logical :: Fortran_Jacobian, box, okgap, oklo, okup, okij, ok_tran, prn, ok
      Character(Len=1) :: flagx, flagt
      Character(Len=10) :: rstr
      Character(Len=20) :: jstr, jfdstr
      Character(Len=200) :: rec
      Continue

      test_tol = iparams%options%derivative_test_tol
      perturbation = iparams%options%fd_step
      Fortran_Jacobian = iparams%options%Fortran_Jacobian
      box = iparams%box%has_box
      prlvl = iparams%options%print_level

      ! Consistency check
      if (Allocated(iparams%f_pert)) Then
         iparams%inform%status = NLLS_ERROR_UNEXPECTED
         goto 100
      end if
      if (.Not. Allocated(iparams%f)) Then
         iparams%inform%status = NLLS_ERROR_UNEXPECTED
         goto 100
      end if

      Allocate(iparams%f_pert(m), stat=ierr)
      if (ierr /= 0) Then
         iparams%inform%status = NLLS_ERROR_ALLOCATION
         goto 100
      end if

      Allocate(J_fd(n * m), Stat=ierr)
      if (ierr /= 0) then
         iparams%inform%status = NLLS_ERROR_ALLOCATION
         goto 100
      end if

      Call fd_jacobian(ierr, n, m, iparams%x, iparams%f, J_fd, iparams)
      if (ierr == -2031 ) then
         ! eval_f provided a rubbish FD estimation
         iparams%inform%external_name = 'eval_F/fd_jacobian'
         iparams%inform%external_return = 2031
         iparams%inform%status = NLLS_ERROR_FROM_EXTERNAL
         goto 100
      elseif (ierr /= 0) then
         ! error from user eval_F callback
         iparams%inform%external_name = 'eval_F'
         iparams%inform%external_return = 2101
         iparams%inform%status = NLLS_ERROR_FROM_EXTERNAL
         goto 100
      end if

      Call printmsg(1,.False.,iparams%options,1, '')
      Write(rec, Fmt=99999)
      Call printmsg(1,.False.,iparams%options,1, rec)
      Call printmsg(1,.False.,iparams%options,1, '')
      Write(rec, Fmt=99997) merge('Fortran (column-major)', &
         'C (row-major)         ', Fortran_Jacobian)
      Call printmsg(1,.False.,iparams%options,1, rec)
      Call printmsg(1,.False.,iparams%options,1, '')

      ierr = 0
      tcnt = 0
      skip = 0
      Do ivar = 1, n
         this_perturbation = perturbation * max(1.0_wp, abs(iparams%x(ivar)));
         okgap = .True.
         oklo = .True.
         okup = .True.
         if (box) then
            okgap = iparams%box%blx(ivar) < iparams%box%bux(ivar)
            oklo = iparams%box%blx(ivar) <= iparams%x(ivar) - this_perturbation
            okup = iparams%x(ivar) + this_perturbation <= iparams%box%bux(ivar)
            okgap = okgap .and. (oklo .or. okup)
         end if

         Do jcon = 1, m
            if (Fortran_Jacobian) then
               idx = (ivar-1) * m + jcon
               idx_tran = (jcon-1) * n + ivar ! get also the transpose
               iivar = ivar ! Fortran indexing
               jjcon = jcon
            else
               idx_tran = (ivar-1) * m + jcon ! get also the transpose
               idx = (jcon-1) * n + ivar
               iivar = ivar - 1 ! C indexing
               jjcon = jcon - 1
            end if
            relerr = abs(J_fd(idx) - J(idx)) / max(abs(J_fd(idx)), test_tol)
            relerr_tran = abs(J_fd(idx) - J(idx_tran)) / max(abs(J_fd(idx)), test_tol)
            okij = relerr <= test_tol
            ok_tran = .False.
            if (okgap) then
               ! Only tally if there is a gap to calculate the FD approximation
               ! and allow to solve a problems with fixed variables
               if (.not. okij) ierr = ierr + 1 ! Tally derivatives to verify/fix
               ok_tran = relerr_tran <= test_tol
               if ((.not. okij) .and. ok_tran) tcnt = tcnt + 1 ! Tally transpose
            else
               skip = skip + 1
            end if
            prn = (.not. okij ) .or. prlvl >= 2
            if (prn) then
               flagx = merge(' ', 'X', okij)
               flagt = merge('T', ' ', (.Not. okij) .And. ok_tran)
               if (J(idx) /= 0.0 .And. abs(J(idx)) < 1.0e-99_wp) then
                  jstr = ' ~0.0                '
               else if (J(idx)> 9.99e+99_wp) then
                  jstr = '       +Inf         '
               else if (J(idx)< -9.99e+99_wp) then
                  jstr = '       -Inf         '
               else
                  Write(jstr, '(Es20.12e2)') J(idx)
               end if
               if (J_fd(idx) /= 0.0 .And. abs(J_fd(idx)) < 1.0e-99_wp) then
                  jfdstr = '~0.0                '
               else if (J_fd(idx)> 9.99e+99_wp) then
                  jfdstr = '       +Inf         '
               else if (J_fd(idx)< -9.99e+99_wp) then
                  jfdstr = '       -Inf         '
               else
                  Write(jfdstr, '(Es20.12e2)') J_fd(idx)
               end if
               if (relerr /= 0.0 .And. relerr < 1.0e-99_wp) then
                  rstr = '~0.0     '
               else if (relerr > 9.9e+99_wp) then
                  rstr = '   +Inf   '
               else
                  write(rstr, '(Es10.3e2)') relerr
               end if
               Write(rec, Fmt=99911) iivar, jjcon, Jstr, Jfdstr, rstr, &
                  this_perturbation, flagx, flagt, merge('    ', 'Skip', okgap)
               Call printmsg(1,.False.,iparams%options,1, rec)
            end if
         End Do
      End Do

      ! Print summary

      Call printmsg(1,.False.,iparams%options,1, '')
      if (ierr == 0) then
         Write(rec, Fmt=80001)
      else
         Write(rec, Fmt=80000) ierr
      end if
      Call printmsg(1,.False.,iparams%options,1, rec)

      if (skip > 0) then
         Call printmsg(1,.False.,iparams%options,1, '')
         Write(rec, Fmt=70000) skip
         Call printmsg(1,.False.,iparams%options,1, rec)
      end if

      if (tcnt > 0) then
         Call printmsg(1,.False.,iparams%options,1, '')
         Write(rec, Fmt=80002) tcnt
         Call printmsg(1,.False.,iparams%options,1, rec)
         Write(rec, Fmt=80003)
         Call printmsg(1,.False.,iparams%options,1, rec)
      end if


      Write(rec, Fmt=99998)
      Call printmsg(1,.False.,iparams%options,1, '')
      Call printmsg(1,.False.,iparams%options,1, rec)
      Call printmsg(1,.False.,iparams%options,1, '')

      if (ierr > 0) iparams%inform%status = NLLS_ERROR_BAD_JACOBIAN

100   continue

      If (Allocated(J_fd)) Deallocate(J_fd)
      If (Allocated(iparams%f_pert)) Deallocate(iparams%f_pert)

99999 Format (1X,'Begin Derivative Checker')
99998 Format (1X,'End Derivative Checker')
99997 Format (4X,'Jacobian storage scheme (Fortran_Jacobian) = ',A)
99911 Format (4X,'Jac[',I6,',',I6,'] = ',A20,' ~ ', A20, 2X,'[',A10,'], (',E10.3e2,')',2X,2(A1),1X,A)
80000 Format (4X,'Derivative checker detected ',I6,1X,'likely error(s)')
80001 Format (4X,'It seems that derivatives are OK.')
80002 Format (4X,'Note: derivative checker detected that ',I6,' entries may correspond to the transpose.')
80003 Format (4X,'Verify the Jacobian storage ordering is correct.')
70000 Format (4X,'Warning: derivative checker skipped ',I6,' entries that have too tight bounds on the variable(s).')

   End Subroutine check_jacobian

   ! FINITE DIFFERENCE DRIVER for the JACOBIAN
   ! This was inspired by IpOpt TNLPAdapter::internal_eval_jac_g
   ! Assumes that eval_F was called previously and that f(:) = f(x), the residual
   ! evaluate at x
   Subroutine fd_jacobian(status, n, m, x, f, J, iparams)
      Implicit None
      integer, intent(out) :: status
      integer, intent(in) :: n, m
      double precision, dimension(*), intent(inout) :: x
      double precision, dimension(*), intent(in) :: f
      double precision, dimension(*), intent(out) :: J
      type(params_internal_type), intent(inout) :: iparams

      integer :: ivar, jcon, idx
      Real(Kind=wp) :: perturbation, this_perturbation, xorig
      Logical :: Fortran_Jacobian, okgap, oklo, okup, box

      Continue

      perturbation = iparams%options%fd_step
      Fortran_Jacobian = iparams%options%Fortran_Jacobian

      box = iparams%box%has_box
      ! Compute the finite difference Jacobian
      ! Forward / Backward step
      Do ivar = 1, n
         this_perturbation = perturbation * max(1.0_wp, abs(x(ivar)));
         okgap = .True.
         oklo = .True.
         okup = .True.
         if (box) then
            okgap = iparams%box%blx(ivar) < iparams%box%bux(ivar)
            oklo = iparams%box%blx(ivar) <= x(ivar) - this_perturbation
            okup = x(ivar) + this_perturbation <= iparams%box%bux(ivar)
            okgap = okgap .and. (oklo .or. okup)
            ! if okgap = FALSE then ith variable is fully constrained, or
            ! bounds are too tight to use...
            ! a sensible step... Set derivatives to zero
         end if

         if (okgap) then
            xorig = x(ivar)
            ! either okup or oklo or both are good to use
            if (.not. okup) this_perturbation = -this_perturbation ! reverse direction
            x(ivar) = x(ivar) + this_perturbation

            Call iparams%user%eval_F(status, n, m, x(1:n), iparams%f_pert(1:m), iparams%user%params)
            iparams%inform%f_eval = iparams%inform%f_eval + 1 ! global F ledger
            iparams%inform%fd_f_eval = iparams%inform%fd_f_eval + 1 ! legder for FD calls
            x(ivar) = xorig
            ! For now no recovery at this point
            if (status /= 0) return
         end if

         Do jcon = 1, m
            if (Fortran_Jacobian) then
               idx = (ivar-1) * m + jcon
            else
               idx = (jcon-1) * n + ivar
            end if
            J(idx) = merge(( iparams%f_pert(jcon) - f(jcon) ) / this_perturbation, 0.0_wp, okgap)
            if (J(idx) /= J(idx)) then
               ! Estimate is rubbish
               status = -2031
               return
            end if
         End Do
      End Do

   End Subroutine fd_jacobian


end module ral_nlls_fd

