Module ral_nlls_bounds
Implicit None
Private
! User routine
Public :: nlls_setup_bounds
! Routines used by ral_nlls_internal
Public :: box_proj, box_projdir

Contains

Subroutine nlls_setup_bounds(params, n, blx, bux, options, inform)
    Use ral_nlls_workspaces
    Implicit None
    Integer, Intent(In)                       :: n
    Type(NLLS_inform), Intent(InOut)          :: inform
    Type(NLLS_options), Intent(In)            :: options
    Class(params_box_type), Intent(InOut)     :: params
    Real(Kind=wp), Intent(In)                 :: blx(n), bux(n)
    Integer                                   :: i, iusrbox, ierr

    Continue

    params%prjchd = .False.
    params%quad_i = 0
    params%quad_c = 0.0_wp
    params%quad_q = 0.0_wp
    
    iusrbox = 0
    Do i = 1, n
      If ( blx(i) <= bux(i) .And. blx(i)==blx(i) .And. bux(i)==bux(i) ) Then
        If ( blx(i) > -options%box_bigbnd .Or. options%box_bigbnd > bux(i)) Then
          iusrbox = 1
        End If
      Else
        inform%status = NLLS_ERROR_BAD_BOX_BOUNDS
        Go To 100
      End If
    End Do
    
    params%iusrbox = iusrbox
    If (iusrbox==1) Then
      params%iusrbox = 1
      ! Clear all arrays...
      If (allocated(params%blx)) deallocate(params%blx, Stat=ierr)
      If (allocated(params%bux)) deallocate(params%bux, Stat=ierr)
      If (allocated(params%pdir)) deallocate(params%pdir, Stat=ierr)
      If (allocated(params%normFref)) deallocate(params%normFref, Stat=ierr)
      If (allocated(params%sk)) deallocate(params%sk, Stat=ierr)
      If (allocated(params%g)) deallocate(params%g, Stat=ierr)
      Allocate(params%blx(n), params%bux(n), params%pdir(n), params%g(n),      &
        params%normFref(options%box_nFref_max), params%sk(n), Stat=ierr)
      if (ierr /= 0) Then
        If (allocated(params%blx)) deallocate(params%blx, Stat=ierr)
        If (allocated(params%bux)) deallocate(params%bux, Stat=ierr)
        If (allocated(params%pdir)) deallocate(params%pdir, Stat=ierr)
        If (allocated(params%normFref)) deallocate(params%normFref, Stat=ierr)
        If (allocated(params%sk)) deallocate(params%sk, Stat=ierr)
        If (allocated(params%g)) deallocate(params%g, Stat=ierr)
        inform%status = NLLS_ERROR_ALLOCATION
        inform%bad_alloc = 'ral_nlls_box'
        Go To 100
      end if
      params%normfref(1:options%box_nFref_max) = -1.0e-20_wp
      Do i = 1, n
        params%blx(i) = max(blx(i), -options%box_bigbnd)
        params%bux(i) = min(bux(i), options%box_bigbnd)
      End Do
    End If

100 Continue

End Subroutine nlls_setup_bounds

Subroutine box_proj(w, n, x, xnew, dir, alpha)
    Use ral_nlls_workspaces, Only: box_type, wp
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
      If (w%iusrbox==1) Then
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

    If (w%iusrbox==1) Then
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
    Use ral_nlls_workspaces, Only: params_box_type, wp, box_type
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
    
    If (w%iusrbox==1) Then
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
  
  
End Module ral_nlls_bounds
