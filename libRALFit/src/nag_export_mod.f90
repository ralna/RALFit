Module NAG_EXPORT_MOD
      Use ral_nlls_workspaces, Only: wp
      Implicit None
      Private
      Public :: e04rlln, e04rlpn, calculate_covm

      Integer, Parameter :: nps_ils_alpnorig = 1
      Integer, Parameter :: nps_ils_alporig = 2
      Integer, Parameter :: nps_ils_bracket = 3
      Integer, Parameter, Public :: nps_ils_it = 4
      Integer, Parameter :: nps_ils_itbad = 5
      Integer, Parameter :: nps_ils_itflat = 6
      Integer, Parameter :: nps_ils_lastbad = 7
      Integer, Parameter :: nps_ils_lastflat = 8
      Integer, Parameter :: nps_ils_maxit = 9
      Integer, Parameter :: nps_ils_status = 10
      Integer, Parameter :: nps_ils_wrnflat = 11
      Integer, Parameter, Public :: nps_ils_len = 12

      Integer, Parameter :: nps_rls_alp = 1
      Integer, Parameter :: nps_rls_alpl = 2
      Integer, Parameter :: nps_rls_alpn = 3
      Integer, Parameter :: nps_rls_alpu = 4
      Integer, Parameter :: nps_rls_dirnrm = 5
      Integer, Parameter :: nps_rls_eta = 6
      Integer, Parameter :: nps_rls_f0 = 7
      Integer, Parameter :: nps_rls_falp = 8
      Integer, Parameter :: nps_rls_fl = 9
      Integer, Parameter :: nps_rls_fu = 10
      Integer, Parameter :: nps_rls_maxalp = 11
      Integer, Parameter :: nps_rls_mu1 = 12
      Integer, Parameter :: nps_rls_mu2 = 13
      Integer, Parameter :: nps_rls_pi0 = 14
      Integer, Parameter, Public :: nps_rls_len = 15

      Real (Kind=wp), Parameter :: x02ajfn_nag = epsilon(1.0E0_wp)/2.0_wp
    Contains

    Subroutine e04rlpn(f0,pi0,dirnrm,alpn,maxalp,eta,mu1,mu2,maxit,ils,rls,    &
      iflag)
!     NAG COPYRIGHT 2015.

!     > e04rlpn/ls_setup - given various constants (options) and
!     > initial conditions, setup data in workspace arrays for various
!     > line search procedures within the suite ls_test/e04rlmn(),
!     > ls_test_simple/e04rlln(), ls_print/e04rlnn().
!     >
!     > Options for LS are stored in ils/rls arrays so that the procedures
!     > are independent of Pennon and can be used elsewhere, however, this
!     > setup routine serves for more than 1 LS so it requires all options
!     > of all LS it serves for. It's not great but...
!     >
!     > Input:
!     > - f0 = f(x+0*d), i.e., line search starting point
!     > - pi0 = grad f(x)*d, i.e., directional derivative
!     > - dirnrm = ||d||, used only for printing (for relative scale of pi0)
!     > - alpn is the next step to try, thus the initial step length
!     >     (first trial step) to be provided in the next call to
!     >     ls_test/e04rlmn(), it should be most likely 1.0 for Newton method;
!     > - maxalp (>=alpn) is the maximal allowed step (no limit if <=0),
!     >     [actually not really used right now as none of the LS extrapolates]
!     > - eta = sufficient decrease constant for Armijo LS
!     > - mu1, mu2 = sufficient decrease and 'not-too-short-step' for Goldstein
!     > - maxit = max steps given by the user but it might be reduced by LS
!     >     if the input seems bad (i.e., fail rather early)
!     >
!     > Output:
!     > - ils/rls - workarrays for all LS
!     > Returns iflag=1 if OK, iflag=-20 if not a descent direction
!     >
!     > \ingroup e04r_solver

!     .. Implicit None Statement ..
      Implicit None
!     .. Scalar Arguments ..
      Real (Kind=wp), Intent (In)      :: alpn, dirnrm, eta, f0, maxalp, mu1,  &
                                          mu2, pi0
      Integer, Intent (Out)            :: iflag
      Integer, Intent (In)             :: maxit
!     .. Array Arguments ..
      Real (Kind=wp), Intent (Out)     :: rls(nps_rls_len)
      Integer, Intent (Out)            :: ils(nps_ils_len)
!     .. Local Scalars ..
      Real (Kind=wp)                   :: epsmach
!     .. Intrinsic Procedures ..
      Intrinsic                        :: abs, max, min
!     .. Executable Statements ..
      Continue

!     OK flag (=requesting evaluation in LS)
      iflag = 1

!     option settings of the line searches
!     typical defaults: eta=0.01, mu1/2 = 0.1, 0.9, maxit=30
      rls(nps_rls_eta) = eta
      rls(nps_rls_mu1) = mu1
      rls(nps_rls_mu2) = mu2
      ils(nps_ils_maxit) = maxit

!     initial conditions
!     max step...
      If (maxalp>0.0_wp) Then
        rls(nps_rls_maxalp) = maxalp
      Else
        rls(nps_rls_maxalp) = 1e10_wp
      End If
!     next step to try (initial step)
      rls(nps_rls_alpn) = alpn
      ils(nps_ils_alpnorig) = 1
!     'current' step (starting point)
      rls(nps_rls_f0) = f0
      rls(nps_rls_alp) = 0.0_wp
      rls(nps_rls_falp) = f0
      ils(nps_ils_alporig) = 0
!     derivative for alpha = 0 (=grad f'*dir)
      rls(nps_rls_pi0) = pi0
      rls(nps_rls_dirnrm) = dirnrm

!     zero counters etc.
      ils(nps_ils_it) = 0
      ils(nps_ils_itbad) = 0
      ils(nps_ils_lastbad) = -1
      ils(nps_ils_itflat) = 0
      ils(nps_ils_lastflat) = -1
      ils(nps_ils_wrnflat) = 0
      ils(nps_ils_status) = 0

!     no bracket available (1 top, 2 bottom, 3 both)
!     TODO once if interval shortened (upperbnd), mark it as top bracket?
      ils(nps_ils_bracket) = 0
      rls(nps_rls_alpu) = 0.0_wp
      rls(nps_rls_fu) = 0.0_wp
      rls(nps_rls_alpl) = 0.0_wp
      rls(nps_rls_fl) = 0.0_wp
!     stopping flags (mainly for printing)
      ils(nps_ils_status) = 0

      If (pi0>=0.0_wp) Then
!       not a descent direction, cannot continue
        iflag = -20
      End If

!     Check if we are close to machine accuracy with function values:
!     Most line searches require a sufficient decrease condition to be
!     satisfied, i.e.,
!     f(alp) <= f(0) + alp*mu*pi0
!     If pi0 is pretty small w.r.t. f0 (either the direction is rubbish
!     or we are very close to the solution), RHS will be practically f0
!     and we would rely on any marginal changes (noise) in f(alp).
!     If 'flat evaluation' is expected, reduce number of LS iterations
!     and let LS now - it can either accept the step based on gradient
!     changes or perform extrapolation or something.
      epsmach = x02ajfn_nag
      If (abs(pi0)<=1E+3_wp*epsmach*abs(f0)) Then
        ils(nps_ils_wrnflat) = 1
        ils(nps_ils_maxit) = min(maxit,max(5,maxit/4))
      End If

      Return
    End Subroutine e04rlpn
    
    Subroutine e04rlln(armijo,evalok,alpn,falp,ils,rls,iflag)
!     NAG COPYRIGHT 2015.

!     > ls_test_simple/e04rlln - main line search procedure for either simple
!     > backtracking Armijo LS or full step acceptance (i.e., shorten
!     > only if you cannot evaluate).
!     >
!     > It is written with RCI (reverse communication interface). Call first
!     > ls_setup/e04rlpn() and suggest the first trial step (e.g., alpn=1.0).
!     > Then call this routine with falp=f(alpn) evaluated at the initial step.
!     > At each call exactly one function evaluation (falp) is processed.
!     > Either the point is rejected and alpn contains on return a new trial
!     > step or if successful the accepted step (also stored in
!     > rls(nps_rls_alp)) and LS ends. If LS failes (e.g., max LS iter),
!     > alpn is the last tried step, i.e., no new trial step is generated.
!     >
!     > See ls_test/e04rlmn() for further details. This is its simple brother.
!     > ls_test/e04rlmn() performs stricter Goldstein tests, interpolation &
!     > safeguarding.
!     >
!     > Input:
!     >   evalok, falp = if (evalok) falp should hold the value of point
!     >     requested last time.
!     >     (intent out of falp is just for convenience, no data is returned)
!     >
!     > workspace ils/rls as setup by ls_setup/e04rlpn(), shouldn't be changed
!     > by user.
!     >
!     > Output:
!     >   alpn = new trial step size or final step size
!     >   iflag=1 ... evaluate falp and call again (at alpn)
!     >   iflag=0 ... finished OK (accept step)
!     >   iflag<0 ... warning (-2 "eval flat" - conditional acceptance)
!     >   iflag<=-10 ... error (-10 max ls it, -99 unexpected)
!     >
!     > \ingroup e04r_solver



!     .. Implicit None Statement ..
      Implicit None
!     .. Scalar Arguments ..
      Real (Kind=wp), Intent (Out)     :: alpn
      Real (Kind=wp), Intent (Inout)   :: falp
      Integer, Intent (Out)            :: iflag
      Logical, Intent (In)             :: armijo, evalok
!     .. Array Arguments ..
      Real (Kind=wp), Intent (Inout)   :: rls(nps_rls_len)
      Integer, Intent (Inout)          :: ils(nps_ils_len)
!     .. Local Scalars ..
      Real (Kind=wp)                   :: alp, epsmach
      Logical                          :: conditional, conv, c_decrease,       &
                                          evalflat
!     .. Intrinsic Procedures ..
      Intrinsic                        :: abs, merge
!     .. Executable Statements ..
      Continue

!     if function cannot be evaluated at the given point,
!     keep backtracking till a good point is found again
      If (.Not. evalok) Then
        falp = 1e10_wp
        ils(nps_ils_itbad) = ils(nps_ils_itbad) + 1
      End If
      iflag = -99

!     new iteration
      ils(nps_ils_it) = ils(nps_ils_it) + 1
!     keep current tested step in alpn in case no new trial
!     is generated (successful stop or hard fail)
      alpn = rls(nps_rls_alpn)
!     alp ... our current trial step, store the value
      alp = rls(nps_rls_alpn)
      rls(nps_rls_alp) = rls(nps_rls_alpn)
      rls(nps_rls_falp) = falp
      ils(nps_ils_alporig) = ils(nps_ils_alpnorig)

!     almost not distinquishable function value from f0?
      epsmach = x02ajfn_nag
      evalflat = abs(falp-rls(nps_rls_f0)) <= 10.0_wp*epsmach*(abs(falp)+abs(  &
        rls(nps_rls_f0)))
      If (evalflat) Then
        ils(nps_ils_itflat) = ils(nps_ils_itflat) + 1
      End If

!     (A) TEST INCOMING POINT

!     sufficient decrease
      c_decrease = falp < rls(nps_rls_f0) + alp*rls(nps_rls_eta)*rls(          &
        nps_rls_pi0)

!     if (armijo), require sufficient decrease condition, otherwise just
!     a successful evaluation and accept it anyway
      conv = evalok .And. (c_decrease .Or. .Not. armijo)

!     conditional acceptance when expected 'flat evaluations' happened
      conditional = evalok .And. ils(nps_ils_wrnflat) /= 0 .And. evalflat

!     set flags for printing - how to convert logicals to bits??
      ils(nps_ils_status) = merge(1,merge(4,2,conditional),conv) +             &
        merge(8,0,c_decrease) + merge(128,0,evalflat) + merge(0,256,evalok)

      If (conv) Then
!       step accepted, stop, (alpn=alp already)
        iflag = 0
        Go To 100
      Else If (conditional) Then
!       conditional acceptance (i.e., test gradients and decide outside of LS)
        iflag = -2
        Go To 100
      Else If (ils(nps_ils_it)>ils(nps_ils_maxit)) Then
        iflag = -10
        Go To 100
      End If

!     (B) GENERATE A TRIAL POINT

      If (.Not. evalok) Then
!       not much to do, just (slowly) shorten the step
        If (ils(nps_ils_itbad)<=3) Then
!         try to do small steps first
          alpn = alp*0.75_wp
        Else
!         then a bit more significant shortening
          alpn = alp*0.5_wp
        End If
        ils(nps_ils_alpnorig) = 6
        iflag = 1
      Else
!       basic Armijo backtracking
        alpn = 0.5_wp*alp
        ils(nps_ils_alpnorig) = 7
        iflag = 1
      End If

!     remember the suggested trial step length in the next call
      rls(nps_rls_alpn) = alpn

100   Continue
      Return
    End Subroutine e04rlln

    Subroutine calculate_covm(m,n,j,inform,options,iflag)
!     Build Covariance Matrix

!     Returns one of the following elements:
!     iflag       Return
!       0         do nothing
!       1         covariance matrix: COV = sigma2 * (J^T J)^-1
!       2         diagonal elements of COV matrix
!       3         product COV = J^T J

!     Notes:
!      * On error cov is not allocated

!     NAG COPYRIGHT 2020.
!     .. Use Statements ..
      Use ral_nlls_workspaces, Only: nlls_options, nlls_inform
!     .. Implicit None Statement ..
      Implicit None
!     .. Scalar Arguments ..
      Type (nlls_options), Intent (In) :: options
      Integer, Intent (In)             :: iflag
      Integer, Intent (In)             :: m, n
!     .. Array Arguments ..
      Type (nlls_inform), Intent (Inout) :: inform
      Real (Kind=wp), Intent (In)      :: j(*)
!     .. Local Scalars ..
      Real (Kind=wp)                   :: sigma2
      Integer                          :: ierr, irank, i, ld
      Character (1)                    :: transa, transb
!     .. Intrinsic Procedures ..
      Intrinsic                        :: allocated, merge, real
!     real(wp) :: dlange
      External dgemm, dpotrf, dpotri
!     .. Executable Statements ..

      Continue

      If (iflag<=0 .Or. inform%obj<=0.0_wp .Or. iflag>3) Then
!       Wrong param inputs, signal cov matrix is not available
!       LCOV_EXCL_START
        If (allocated(inform%cov)) Then
          Deallocate(inform%cov, stat=ierr)
        End If
        Go To 100
!       LCOV_EXCL_STOP
      End If

      Allocate(inform%cov(n, n), stat=ierr)
      If (ierr/=0) Then
!       LCOV_EXCL_START
        Go To 100
!       LCOV_EXCL_STOP
      End If

!     inform%cov(:,:) = 0.0_wp
      If (options%fortran_jacobian) Then
        transa = 'T'
        transb = 'N'
        ld = m
      Else
        transa = 'N'
        transb = 'T'
        ld = n
      End if

!     Build cov = J^T J
      Call dgemm(transa, transb, n, n, m, 1.0_wp, j, ld, j, ld, 0.0_wp, inform%cov, n)
      If (iflag==3) Then
!       Return  only cov = H = J^T J
        Go To 100
      End If

!     Invert J^T J using LAPACK
      Call dpotrf('u', n, inform%cov, n, ierr)
      if (ierr /= 0) Then
!       Cholesky factorization failed
        Deallocate(inform%cov, stat=ierr)
        Go To 100
      End If

      Call dpotri('u', n, inform%cov, n, ierr)
      if (ierr /= 0) Then
!       Inverse of Cholesky factor failed
        Deallocate(inform%cov, stat=ierr)
        Go To 100
      End If

!     If Cholesky does not breakdown, then
      irank = n
      If (irank < m) Then
        sigma2 = 2.0_wp*inform%obj/real(m-irank,kind=wp)
!       multiply by sigma^2 = 2*obj / (m-k)
        inform%cov(:,:) = sigma2*inform%cov(:,:)
      Else
!       covariance not available for underdetermined system
        Deallocate(inform%cov, stat=ierr)
        Go To 100
      End If

      If (iflag==2) Then
        If (Allocated(inform%var)) Then
          Deallocate(inform%var, stat=ierr)
          If (ierr/=0) Then
!           LCOV_EXCL_START
!           Allocation error
            Deallocate(inform%cov, stat=ierr)
            Go To 100
!           LCOV_EXCL_STOP
          End If
        End If
        Allocate(inform%var(n), stat=ierr)
        If (ierr/=0) Then
!         LCOV_EXCL_START
!         Allocation error
          Deallocate(inform%cov, stat=ierr)
          Go To 100
!         LCOV_EXCL_STOP
        End If
!       Copy only diagonal elements
        Do i = 1, n
          inform%var(i) = inform%cov(i,i)
        End Do 
        Deallocate(inform%cov, stat=ierr)
      End If

100   Continue

    End Subroutine calculate_covm

End Module NAG_EXPORT_MOD
