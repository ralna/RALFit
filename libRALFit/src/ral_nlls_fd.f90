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

   Public :: eval_f_wrap, eval_j_wrap, eval_hf_wrap, eval_hp_wrap
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
   Type, Public :: jacobian_handle
      Logical :: setup = .False.
      Integer :: n, m
      Procedure(eval_f_type), NoPass, Pointer :: eval_f
      Class(params_base_type), Pointer :: params
      Real(Kind=wp), Pointer :: lower, upper
      Logical :: f_storage
      ! private elements for internal params
      Type(params_internal_type) :: iparams
      Type(nlls_inform) :: inform
      Type(NLLS_options) :: options
      Type(box_type) :: box
      Real(Kind=wp), Dimension(:), Allocatable :: f
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
      Use ral_nlls_workspaces, Only : params_internal_type, wp, NLLS_options, nlls_inform, setup_bounds_type, setup_iparams_type
      Implicit None
      Integer, Intent(out) :: status
      Type(jacobian_handle), Intent(Inout) :: handle
      Integer, Intent(in) :: n, m
      Real(Kind=wp), Dimension(n) :: x
      Procedure( eval_f_type ), Pointer :: eval_f
      Class(params_base_type), Intent(inout), Optional :: params
      Real(Kind=wp), Dimension(n), Intent(in), Target, Optional :: lower, upper
      Logical, Optional :: f_storage

      Type(params_base_type) :: no_params
      Integer :: ierr_dummy

      Continue

      Allocate(handle%f(m), Stat=status)
      If (status /= 0) Then
         status = -99
         GoTo 100
      End If

      Allocate(handle%iparams%f_pert(m), Stat=status)
      If (status /= 0) Then
         status = -99
         GoTo 100
      End If

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

      If (Present(params)) Then
         Call setup_iparams_type(handle%iparams, params, eval_f,   &
            ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, handle%inform, handle%options, handle%box, x)
      Else
         Call setup_iparams_type(handle%iparams, no_params, eval_f,   &
            ral_nlls_eval_j_dummy, ral_nlls_eval_hf_dummy, handle%inform, handle%options, handle%box, x)
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
      If (Allocated(handle%f)) Deallocate(handle%f)
      Call free_iparams_type(handle%iparams)
      Call remove_workspace_bounds(handle%box)
   End Subroutine jacobian_free

   subroutine jacobian_calc(status, handle, x, J, fd_step)
      use ral_nlls_workspaces, only : wp
      implicit none
      integer, intent(out) :: status
      Type(jacobian_handle), Intent(Inout) :: handle
      Real(Kind=wp), dimension(handle%n), intent(in)  :: x
      Real(Kind=wp), dimension(handle%m*handle%n), intent(out) :: J
      Real(Kind=wp), optional :: fd_step

      Continue

      If (.Not. handle%setup) Then
         status =  -199
         return
      End If

      If (present(fd_step)) Then
         handle%iparams%options%fd_step = fd_step
      End If

      Call fd_jacobian(status, handle%n, handle%m, J, handle%iparams)

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
               ! User provided a custom eval_J return user status and J values
               return
            end if
         end if
         if ( params%fd_type /= 'N' ) then
            ! Finite-differences
            ! Assumes that eval_F was called previously and that f(:) is valid
            Call fd_jacobian(status, n, m, J, params)
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


   ! FINITE DIFFERENCE DRIVER for the JACOBIAN
   ! This is a Fortran version of IpOpt TNLPAdapter::internal_eval_jac_g
   Subroutine fd_jacobian(status, n, m, J, iparams)
      Implicit None
      integer, intent(out) :: status
      integer, intent(in) :: n, m
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
         this_perturbation = perturbation * max(1.0_wp, abs(iparams%x(ivar)));
         okgap = .True.
         oklo = .True.
         okup = .True.
         if (box) then
            okgap = iparams%box%blx(ivar) < iparams%box%bux(ivar)
            oklo = iparams%box%blx(ivar) <= iparams%x(ivar) - this_perturbation
            okup = iparams%x(ivar) + this_perturbation <= iparams%box%bux(ivar)
            okgap = okgap .and. (oklo .or. okup)
            ! if okgap = FALSE then ith variable is fully constrained, or bounds a too tight to use
            ! a sensible step... Set derivatives to zero
         end if

         if (okgap) then
            xorig = iparams%x(ivar)
            ! either okup or oklo or both are good to use
            iparams%x(ivar) = iparams%x(ivar) + merge(this_perturbation, -this_perturbation, okup)

            Call iparams%user%eval_F(status, n, m, iparams%x, iparams%f_pert, iparams%user%params)
            iparams%inform%f_eval = iparams%inform%f_eval + 1 ! gloval F ledger
            iparams%inform%fd_f_eval = iparams%inform%fd_f_eval + 1 ! legder for FD calls
            iparams%x(ivar) = xorig
            ! For now no recovery at this point
            if (status /= 0) return
         end if

         Do jcon = 1, m
            if (Fortran_Jacobian) then
               idx = (ivar-1) * m + jcon
            else
               idx = (jcon-1) * n + ivar
            end if
            J(idx) = merge(( iparams%f_pert(jcon) - iparams%f(jcon) ) / this_perturbation, 0.0_wp, okgap)
            if (J(idx) /= J(idx)) then
               ! Estimate is rubbish
               status = -2031
               return
            end if
         End Do
      End Do

   End Subroutine fd_jacobian


end module ral_nlls_fd

