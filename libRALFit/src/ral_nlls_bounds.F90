! Copyright (c) 2020, The Numerical Algorithms Group Ltd (NAG)
! All rights reserved.
! Copyright (c) 2020, The Science and Technology Facilities Council (STFC)
! All rights reserved.

#include "preprocessor.FPP"

Module MODULE_PREC(ral_nlls_bounds)
Implicit None
Private
! Routines used by ral_nlls_internal
Public :: box_proj, box_projdir

Contains

Subroutine box_proj(w, n, x, xnew, dir, alpha)
    Use MODULE_PREC(ral_nlls_workspaces), Only: box_type, wp
!   Two modes
!   If xnew and d are present, then project x+alpha*dir, otherwise just
!   make x feasible (ignoring either dir or xnew) and return
!   In either case flag in params%prjchd if the projection altered any entry
    Implicit None
    Integer, Intent(In)                    :: n
    type(box_type), Intent(InOut)          :: w
    Real(Kind=wp), Intent(InOut)           :: x(n)
    Real(Kind=wp), Intent(InOut), Optional :: xnew(n)
    Real(Kind=wp), Intent(In), Optional    :: dir(n), alpha
    Integer                                :: i
    Real(Kind=wp)                          :: alp, xi

    Continue
    w%prjchd = .False.
    If (.not. (present(xnew) .And. present(dir))) Then
!     Make feasible point x and return
      If (w%has_box) Then
        Do i = 1, n
          xi = x(i)
          x(i) = max(min(w%bux(i), x(i)), w%blx(i))
          If (xi/=x(i)) Then
            w%prjchd = .True.
          End If
        End Do
      End If
      Go To 100
    End If

    If (present(alpha)) Then
      alp = alpha
    Else
      alp = 1.0_wp
    End If

    If (w%has_box) Then
      Do i = 1, n
        xnew(i) = max(min(w%bux(i), x(i)+alp*dir(i)), w%blx(i))
        If (xnew(i) /= x(i)+dir(i)) Then
          w%prjchd = .True.
        End If
      End Do
    Else
      xnew(1:n) = x(1:n) + alp* dir(1:n)
    End If

100 Continue

  End Subroutine box_proj

  Subroutine box_projdir(w, n, x, dir, normg, sigma)
    Use MODULE_PREC(ral_nlls_workspaces), Only: wp, box_type
    !   Calculate the projected dir and it's two-norm
    !   Assumes dir = -fdx
    !   If there is no box, then normPD=normg and if pdir is allocated then
    !   copy dir to it.
    Implicit None
    type( box_type ), Intent(InOut)       :: w
    Integer, Intent(In)                   :: n
    Real(Kind=wp), Intent(In)             :: x(n), dir(n), normg
    Real(Kind=wp), Intent(In), Optional   :: sigma

    Real(Kind=wp)                         :: alpb
    Integer                               :: i

    If (w%has_box) Then
       If (Present(sigma)) Then
          alpb = sigma
       Else
          alpb = 1.0_wp
       End If
       w%normPD = 0.0_wp
       do i = 1, n
          If (w%bux(i)/=w%blx(i)) Then
             w%pdir(i) = max(min(w%bux(i), x(i)+alpb*dir(i)), w%blx(i))-x(i)
             w%normPD = w%normPD + (w%pdir(i))**2
          Else
             w%pdir(i) = 0.0_wp
          End If
       end do
       w%normPD = sqrt(w%normPD)
    Else
       If (allocated(w%pdir)) Then
          w%pdir(1:n) = dir(1:n)
       End If
       w%normPD = normg
    End If
  End Subroutine box_projdir


End Module MODULE_PREC(ral_nlls_bounds)
