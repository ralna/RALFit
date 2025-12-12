! Contains :: ral_nlls_roots
!             ral_dtrs
!             ral_drqs
! THIS VERSION: RAL_NLLS 1.1 - 07/03/2016 AT 09:45 GMT.

!-*-*-*-*-*-*-*-*-*-  R A L _ N L L S _ R O O T S   M O D U L E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for RAL_NLLS productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package ROOTS, March 7th, 2016

!  For full documentation, see
!   http://galahad.rl.ac.uk/galahad-www/specs.html

#include "preprocessor.FPP"

   MODULE MODULE_PREC(RAL_NLLS_ROOTS)

!     --------------------------------------------------------------------
!     |                                                                  |
!     |  Find (all the) real roots of polynomials with real coefficients |
!     |                                                                  |
!     --------------------------------------------------------------------

      USE MODULE_PREC(RAL_NLLS_TYPES), ONLY: lp, np, wp
      USE MODULE_PREC(RAL_NLLS_SYMBOLS)

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: ROOTS_quadratic, ROOTS_cubic, ROOTS_quartic

!--------------------
!   P r e c i s i o n
!--------------------

!     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: out = 6
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: four = 4.0_wp
      REAL ( KIND = wp ), PARAMETER :: six = 6.0_wp
      REAL ( KIND = wp ), PARAMETER :: quarter = 0.25_wp
      REAL ( KIND = wp ), PARAMETER :: threequarters = 0.75_wp
      REAL ( KIND = wp ), PARAMETER :: onesixth = one / six
      REAL ( KIND = wp ), PARAMETER :: onethird = one / three
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: twothirds = two / three
!     REAL ( KIND = wp ), PARAMETER :: pi = four * ATAN( 1.0_wp )
      REAL ( KIND = wp ), PARAMETER :: pi = 3.1415926535897931_wp
!     REAL ( KIND = wp ), PARAMETER :: magic = twothirds * pi
      REAL ( KIND = wp ), PARAMETER :: magic = 2.0943951023931953_wp  !! 2 pi/3
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one )

!  interface to LAPACK: eigenvalues of a Hessenberg matrix

      INTERFACE HSEQR
        SUBROUTINE SHSEQR( job, compz, n, ilo, ihi, H, ldh,  WR, WI, Z, ldz,   &
                           WORK, lwork, info )
        import :: lp
        INTEGER, INTENT( IN ) :: ihi, ilo, ldh, ldz, lwork, n
        INTEGER, INTENT( OUT ) :: info
        CHARACTER ( LEN = 1 ), INTENT( IN ) :: compz, job
        REAL(KIND=lp), INTENT( INOUT ) :: H( ldh, * ), Z( ldz, * )
        REAL(KIND=lp), INTENT( OUT ) :: WI( * ), WR( * ), WORK( * )
        END SUBROUTINE SHSEQR

        SUBROUTINE DHSEQR( job, compz, n, ilo, ihi, H, ldh,  WR, WI, Z, ldz,   &
                           WORK, lwork, info )
        import :: np
        INTEGER, INTENT( IN ) :: ihi, ilo, ldh, ldz, lwork, n
        INTEGER, INTENT( OUT ) :: info
        CHARACTER ( LEN = 1 ), INTENT( IN ) :: compz, job
        REAL(KIND=np), INTENT( INOUT ) :: H( ldh, * ), Z( ldz, * )
        REAL(KIND=np), INTENT( OUT ) :: WI( * ), WR( * ), WORK( * )
        END SUBROUTINE DHSEQR
      END INTERFACE HSEQR

   CONTAINS

!-*-*-*-*-*-   R O O T S _ q u a d r a t i c  S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2, debug )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the quadratic equation
!
!                   a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1 and a2 are real
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a2, a1, a0, tol
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2
      LOGICAL, INTENT( IN ) :: debug

!  Local variables

      REAL ( KIND = wp ) :: rhs, d, p, pprime

      rhs = tol * a1 * a1
      IF ( ABS( a0 * a2 ) > rhs ) THEN  !  really is quadratic
        root2 = a1 * a1 - four * a2 * a0
        IF ( ABS( root2 ) <= ( epsmch * a1 ) ** 2 ) THEN ! numerical double root
          nroots = 2 ; root1 = -  half * a1 / a2 ; root2 = root1
        ELSE IF ( root2 < zero ) THEN    ! complex not real roots
          nroots = 0 ; root1 = zero ; root2 = zero
        ELSE                             ! distint real roots
          d = - half * ( a1 + SIGN( SQRT( root2 ), a1 ) )
          nroots = 2 ; root1 = d / a2 ; root2 = a0 / d
          IF ( root1 > root2 ) THEN
            d = root1 ; root1 = root2 ; root2 = d
          END IF
        END IF
      ELSE IF ( a2 == zero ) THEN
        IF ( a1 == zero ) THEN
          IF ( a0 == zero ) THEN         ! the function is zero
            nroots = 1 ; root1 = zero ; root2 = zero
          ELSE                           ! the function is constant
            nroots = 0 ; root1 = zero ; root2 = zero
          END IF
        ELSE                             ! the function is linear
          nroots = 1 ; root1 = - a0 / a1 ; root2 = zero
        END IF
      ELSE                               ! very ill-conditioned quadratic
        nroots = 2
        IF ( - a1 / a2 > zero ) THEN
          root1 = zero ; root2 = - a1 / a2
        ELSE
          root1 = - a1 / a2 ; root2 = zero
        END IF
      END IF

!  perfom a Newton iteration to ensure that the roots are accurate

      IF ( nroots >= 1 ) THEN
        p = ( a2 * root1 + a1 ) * root1 + a0
        pprime = two * a2 * root1 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime
          root1 = root1 - p / pprime
          p = ( a2 * root1 + a1 ) * root1 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 1, root1, p
        IF ( nroots == 2 ) THEN
          p = ( a2 * root2 + a1 ) * root2 + a0
          pprime = two * a2 * root2 + a1
          IF ( pprime /= zero ) THEN
            IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime
            root2 = root2 - p / pprime
            p = ( a2 * root2 + a1 ) * root2 + a0
          END IF
          IF ( debug ) WRITE( out, 2010 ) 2, root2, p
        END IF
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quadratic = ', ES12.4,     &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quadratic = ', ES12.4 )


!  End of subroutine ROOTS_quadratic

      END SUBROUTINE ROOTS_quadratic

!-*-*-*-*-*-*-*-   R O O T S _ c u b i c  S U B R O U T I N E   -*-*-*-*-*-*-*-

      SUBROUTINE ROOTS_cubic( a0, a1, a2, a3, tol, nroots, root1, root2,       &
                              root3, debug )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the cubicc equation
!
!                a3 * x**3 + a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1, a2 and a3 are real
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a3, a2, a1, a0, tol
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2, root3
      LOGICAL, INTENT( IN ) :: debug

!  Local variables

      INTEGER :: info, nroots_q
      REAL ( KIND = wp ) :: a, b, c, d, e, f, p, q, s, t, w, x, y, z
      REAL ( KIND = wp ) :: c0, c1, c2, b0, b1, pprime, u1, u2
      REAL ( KIND = wp ) :: H( 3, 3 ), ER( 3 ), EI( 3 ), ZZ( 1, 3 ), WORK( 33 )

!  define method used:
!    1 = Nonweiler, 2 = Littlewood, 3 = Viete, other = companion matrix

      INTEGER, PARAMETER :: method = 1

!  Check to see if the quartic is actually a cubic

      IF ( a3 == zero ) THEN
        CALL ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2, debug )
        root3 = infinity
        RETURN
      END IF

!  Deflate the polnomial if the trailing coefficient is zero

      IF ( a0 == zero ) THEN
        root1 = zero
        CALL ROOTS_quadratic( a1, a2, a3, tol, nroots, root2, root3, debug )
        nroots = nroots + 1
        RETURN
      END IF

!  1. Use Nonweiler's method (CACM 11:4, 1968, pp269)

      IF ( method == 1 ) THEN
        c0 = a0 / a3
        c1 = a1 / a3
        c2 = a2 / a3

        s = c2 / three
        t = s * c2
        b = 0.5_wp * ( s * ( twothirds * t - c1 ) + c0 )
        t = ( t - c1 ) / three
        c = t * t * t ; d = b * b - c

! 1 real + 2 equal real or 2 complex roots

        IF ( d >= zero ) THEN
          d = ( SQRT( d ) + ABS( b ) ) ** onethird
          IF ( d /= zero ) then
            IF ( b > zero ) then
              b = - d
            ELSE
              b = d
            END IF
            c = t / b
          END IF
          d = SQRT( threequarters ) * ( b - c )
          b = b + c ; c = - 0.5 * b - s
          root1 = b - s
          IF ( d == zero ) THEN
            nroots = 3 ; root2 = c ; root3 = c
          ELSE
            nroots = 1
          END IF

! 3 real roots

        ELSE
          IF ( b == zero ) THEN
            d = twothirds * ATAN( one )
          ELSE
            d = ATAN( SQRT( - d ) / ABS( b ) ) / three
          END IF
          IF ( b < zero ) THEN
            b = two * SQRT( t )
          ELSE
            b = - two * SQRT( t )
          END IF
          c = COS( d ) * b
          t = - SQRT( threequarters ) * SIN( d ) * b - half * c
          d = - t - c - s ; c = c - s ; t = t - s
          IF ( ABS( c ) > ABS( t ) ) then
            root3 = c
          ELSE
            root3 = t
            t = c
          END IF
          IF ( ABS( d ) > ABS( t ) ) THEN
            root2 = d
          ELSE
            root2 = t
            t = d
          END IF
          root1 = t ; nroots = 3
        END IF

!  2. Use Littlewood's method

      ELSE IF ( method == 2 ) THEN
        c2 = a2 / ( three * a3 ) ; c1 = a1 / ( three * a3 ) ; c0 = a0 / a3
        x = c1 - c2 * c2
        y = c0 - c2* ( x + x + c1 )
        z = y ** 2 + four * x ** 3

!  there are three real roots

        IF ( z < zero ) THEN
          a = - two * SQRT( - x )
          b = y / ( a * x )
          y = ATAN2( SQRT( one - b ), SQRT( one + b ) ) * twothirds
          IF ( c2 < zero ) y = y + magic

!  calculate root which does not involve cancellation

          nroots = 1 ; root1 = a * COS( y ) - c2

!  there may be only one real root

        ELSE
          a = SQRT( z ) ; b = half * ( ABS( y ) + a ) ; c = b ** onethird
          IF ( c <= zero ) THEN
            nroots = 3 ; root1 = - c2 ; root2 = - c2 ; root3 = - c2
            GO TO 900
          ELSE
            nroots = 1
            c = c - ( c ** 3 - b ) / ( three * c * c )
            e = c * c + ABS( x )
            f = one / ( ( x / c ) ** 2 + e )
            IF ( x >= zero ) THEN
              x = e / c ; z = y * f
            ELSE
              x = a * f ; z = SIGN( one, y ) * e / c
            END IF
            IF ( z * c2 >= zero ) THEN
              root1 = - z - c2
            ELSE
              root2 = half * z - c2
              root3 = half * SQRT( three ) * ABS( x )
              root1 = - c0 / ( root2 * root2 + root3 * root3 )
              GO TO 900
            END IF
          END IF
        END IF

!  deflate cubic

        b0 = - c0 / root1
        IF ( ABS( root1 ** 3 ) <= ABS( c0 ) ) THEN
          b1 = root1 + three * c2
        ELSE
          b1 = ( b0 - three * c1 ) / root1
        END IF
        CALL ROOTS_quadratic( b0, b1, one, epsmch, nroots_q,                   &
                              root2, root3, debug )
        nroots = nroots + nroots_q


!  3. Use Viete's method

      ELSE IF ( method == 3 ) THEN
        w = a2 / ( three * a3 )
        p = ( a1 / ( three * a3 ) - w ** 2 ) ** 3
        q = - half * ( two * w ** 3 - ( a1 * w - a0 ) / a3 )
        d = p + q ** 2

!  three real roots

        IF ( d < zero ) THEN
          s = ACOS( MIN( one, MAX( - one, q / SQRT( - p ) ) ) )
          p = two * ( - p ) ** onesixth
          nroots = 3
          root1 = p * COS( onethird * ( s + two * pi ) ) - w
          root2 = p * COS( onethird * ( s + four * pi ) ) - w
          root3 = p * COS( onethird * ( s + six * pi ) ) - w

!  one real root

        ELSE
          d = SQRT( d ) ; u1 = q + d ; u2 = q - d
          nroots = 1
          root1 = SIGN( ABS( u1 ) ** onethird, u1 ) +                          &
                  SIGN( ABS( u2 ) ** onethird, u2 ) - w
        END IF

!  4. Compute the roots as the eigenvalues of the relevant compainion matrix

      ELSE
        H( 1, 1 ) = zero ; H( 2, 1 ) = one ; H( 3, 1 ) = zero
        H( 1, 2 ) = zero ; H( 2, 2 ) = zero ; H( 3, 2 ) = one
        H( 1, 3 ) = - a0 / a3 ; H( 2, 3 ) = - a1 / a3 ; H( 3, 3 ) = - a2 / a3
        CALL HSEQR( 'E', 'N', 3, 1, 3, H, 3, ER, EI, ZZ, 1, WORK, 33, info )
        IF ( info /= 0 ) THEN
          IF ( debug ) WRITE( out,                                             &
         &   "( ' ** error return ', I0, ' from HSEQR in ROOTS_cubic' )" ) info
          nroots = 0
          RETURN
        END IF

!  count and record the roots

        nroots = COUNT( ABS( EI ) <= epsmch )
        IF ( nroots == 1 ) THEN
          IF (  ABS( EI( 1 ) ) <= epsmch ) THEN
            root1 = ER( 1 )
          ELSE IF (  ABS( EI( 2 ) ) <= epsmch ) THEN
            root1 = ER( 2 )
          ELSE
            root1 = ER( 3 )
          END IF
        ELSE
          root1 = ER( 1 ) ;  root2 = ER( 2 ) ;  root3 = ER( 3 )
        END IF
      END IF

!  reorder the roots

  900 CONTINUE
      IF ( nroots == 3 ) THEN
        IF ( root1 > root2 ) THEN
          a = root2 ; root2 = root1 ; root1 = a
        END IF
        IF ( root2 > root3 ) THEN
          a = root3
          IF ( root1 > root3 ) THEN
            a = root1 ; root1 = root3
          END IF
          root3 = root2 ; root2 = a
        END IF
        IF ( debug ) WRITE( out, "( ' 3 real roots ' )" )
      ELSE IF ( nroots == 2 ) THEN
        IF ( debug ) WRITE( out, "( ' 2 real roots ' )" )
      ELSE
        IF ( debug ) WRITE( out, "( ' 1 real root ' )" )
      END IF

!  perfom a Newton iteration to ensure that the roots are accurate

      p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      pprime = ( three * a3 * root1 + two * a2 ) * root1 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime
        root1 = root1 - p / pprime
        p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 1, root1, p

      IF ( nroots == 3 ) THEN
        p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        pprime = ( three * a3 * root2 + two * a2 ) * root2 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime
          root2 = root2 - p / pprime
          p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 2, root2, p

        p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        pprime = ( three * a3 * root3 + two * a2 ) * root3 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 3, root3, p, - p / pprime
          root3 = root3 - p / pprime
          p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 3, root3, p
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4,         &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4 )


!  End of subroutine ROOTS_cubic

      END SUBROUTINE ROOTS_cubic

!-*-*-*-*-*-*-   R O O T S _ q u a r t i c   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE ROOTS_quartic( a0, a1, a2, a3, a4, tol, nroots, root1, root2, &
                                root3, root4, debug )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the quartic equation
!
!        a4 * x**4 + a3 * x**3 + a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1, a2, a3 and a4 are real
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a4, a3, a2, a1, a0, tol
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2, root3, root4
      LOGICAL, INTENT( IN ) :: debug

!  Local variables

      INTEGER :: type_roots, nrootsc
      REAL ( KIND = wp ) :: a, alpha, b, beta, c, d, delta, gamma, r
      REAL ( KIND = wp ) :: x1, xm, xmd, xn, xnd
      REAL ( KIND = wp ) :: d3, d2, d1, d0, b4, b3, b2, b1
      REAL ( KIND = wp ) :: rootc1, rootc2, rootc3, p, pprime

!  Check to see if the quartic is actually a cubic

      IF ( a4 == zero ) THEN
        CALL ROOTS_cubic( a0, a1, a2, a3, tol, nroots, root1, root2, root3,    &
                          debug )
        root4 = infinity
        RETURN
      END IF

!  Use Ferrari's algorithm

!  Initialize

      nroots = 0
      b1 = a3 / a4
      b2 = a2 / a4
      b3 = a1 / a4
      b4 = a0 / a4
      d3 = one
      d2 =  - b2
      d1 = b1 * b3 - four * b4
      d0 = b4 * ( four * b2 - b1 * b1 ) - b3 * b3

!  Compute the roots of the auxiliary cubic

      CALL ROOTS_cubic( d0, d1, d2, d3, tol, nrootsc, rootc1, rootc2, rootc3, &
                        debug )
      IF ( nrootsc > 1 ) rootc1 = rootc3
      x1 = b1 * b1 * quarter - b2 + rootc1
      IF ( x1 < zero ) THEN
        xmd = SQRT( - x1 )
        xnd = quarter * ( two * b3 - b1 * rootc1 ) / xmd
        alpha = half * b1 * b1 - rootc1 - b2
        beta = four * xnd - b1 * xmd
        r = SQRT( alpha * alpha + beta * beta )
        gamma = SQRT( half * ( alpha + r ) )
        IF ( gamma == zero ) THEN
          delta = SQRT( - alpha )
        ELSE
          delta = beta * half / gamma
        END IF
        root1 = half * ( - half * b1 + gamma )
        root2 = half * ( xmd + delta )
        root3 = half * ( - half * b1 - gamma )
        root4 = half * ( xmd - delta )
        GO TO 900
      END IF
      IF ( x1 /= zero ) THEN
        xm = SQRT( x1 )
        xn = quarter * ( b1 * rootc1 - two * b3 ) / xm
      ELSE
        xm = zero
        xn = SQRT( quarter * rootc1 * rootc1 - b4 )
      END IF
      alpha = half * b1 * b1 - rootc1 - b2
      beta = four * xn - b1 * xm
      gamma = alpha + beta
      delta = alpha - beta
      a = - half * b1

!  Compute how many real roots there are

      type_roots = 1
      IF ( gamma >= zero ) THEN
        nroots = nroots + 2
        type_roots = 0
        gamma = SQRT( gamma )
      ELSE
        gamma = SQRT( - gamma )
      END IF
      IF ( delta >= zero ) THEN
        nroots = nroots + 2
        delta = SQRT( delta )
      ELSE
        delta = SQRT( - delta )
      END IF
      type_roots = nroots + type_roots

!  Two real roots

      IF ( type_roots == 3 ) THEN
        root1 = half * ( a - xm - delta )
        root2 = half * ( a - xm + delta )
        root3 = half * ( a + xm )
        root4 = half * gamma
        GO TO 900
      ELSE IF ( type_roots /= 4 ) THEN
        IF ( type_roots == 2 ) THEN
          root1 = half * ( a + xm - gamma )
          root2 = half * ( a + xm + gamma )
        ELSE

!  No real roots

          root1 = half * ( a + xm )
          root2 = half * gamma
        END IF
        root3 = half * ( a - xm ) * half
        root4 = half * delta
        GO TO 900
      END IF

!  Four real roots

      b = half * ( a + xm + gamma )
      d = half * ( a - xm + delta )
      c = half * ( a - xm - delta )
      a = half * ( a + xm - gamma )

!  Sort the roots

      root1 = MIN( a, b, c, d )
      root4 = MAX( a, b, c, d )

      IF ( a == root1 ) THEN
        root2 = MIN( b, c, d )
      ELSE IF ( b == root1 ) THEN
        root2 = MIN( a, c, d )
      ELSE IF ( c == root1 ) THEN
        root2 = MIN( a, b, d )
      ELSE
        root2 = MIN( a, b, c )
      END IF

      IF ( a == root4 ) THEN
        root3 = MAX( b, c, d )
      ELSE IF ( b == root4 ) THEN
        root3 = MAX( a, c, d )
      ELSE IF ( c == root4 ) THEN
        root3 = MAX( a, b, d )
      ELSE
        root3 = MAX( a, b, c )
      END IF

  900 CONTINUE

!  Perfom a Newton iteration to ensure that the roots are accurate

      IF ( debug ) THEN
        IF ( nroots == 0 ) THEN
          WRITE( out, "( ' no real roots ' )" )
        ELSE IF ( nroots == 2 ) THEN
          WRITE( out, "( ' 2 real roots ' )" )
        ELSE IF ( nroots == 4 ) THEN
          WRITE( out, "( ' 4 real roots ' )" )
        END IF
      END IF
      IF ( nroots == 0 ) RETURN

      p = ( ( ( a4 * root1 + a3 ) * root1 + a2 ) * root1 + a1 ) * root1 + a0
      pprime = ( ( four * a4 * root1 + three * a3 ) * root1 + two * a2 )       &
                 * root1 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime
        root1 = root1 - p / pprime
        p = ( ( ( a4 * root1 + a3 ) * root1 + a2 ) * root1 + a1 ) * root1 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 1, root1, p

      p = ( ( ( a4 * root2 + a3 ) * root2 + a2 ) * root2 + a1 ) * root2 + a0
      pprime = ( ( four * a4 * root2 + three * a3 ) * root2 + two * a2 )       &
                 * root2 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime
        root2 = root2 - p / pprime
        p = ( ( ( a4 * root2 + a3 ) * root2 + a2 ) * root2 + a1 ) * root2 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 2, root2, p

      IF ( nroots == 4 ) THEN
        p = ( ( ( a4 * root3 + a3 ) * root3 + a2 ) * root3 + a1 ) * root3 + a0
        pprime = ( ( four * a4 * root3 + three * a3 ) * root3 + two * a2 )     &
                   * root3 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 3, root3, p, - p / pprime
          root3 = root3 - p / pprime
          p = ( ( ( a4 * root3 + a3 ) * root3 + a2 ) * root3 + a1 ) * root3 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 3, root3, p

        p = ( ( ( a4 * root4 + a3 ) * root4 + a2 ) * root4 + a1 ) * root4 + a0
        pprime = ( ( four * a4 * root4 + three * a3 ) * root4 + two * a2 )     &
                   * root4 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 4, root4, p, - p / pprime
          root4 = root4 - p / pprime
          p = ( ( ( a4 * root4 + a3 ) * root4 + a2 ) * root4 + a1 ) * root4 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 4, root4, p
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quartic = ', ES12.4,       &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quartic = ', ES12.4 )

!  End of subroutine ROOTS_quartic

      END SUBROUTINE ROOTS_quartic

!  End of module ROOTS

   END MODULE MODULE_PREC(RAL_NLLS_ROOTS)

! THIS VERSION: RAL_NLLS 1.2 - 04/03/2020 AT 13:15 GMT.

!-*-*-*-*-*-*-*-  R A L _ N L L S _ D T R S  double  M O D U L E  *-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package TRS, December 22nd, 2015

   MODULE MODULE_PREC(RAL_NLLS_DTRS)

!       -----------------------------------------------
!      |                                               |
!      | Solve the trust-region subproblem             |
!      |                                               |
!      |    minimize     1/2 <x, H x> + <c, x> + f     |
!      |    subject to      ||x||_2 <= radius          |
!      |    or              ||x||_2  = radius          |
!      |                                               |
!      | where H is diagonal                           |
!      |                                               |
!       -----------------------------------------------

      USE MODULE_PREC(ral_nlls_types), Only: lp, np, wp
      USE MODULE_PREC(RAL_NLLS_SYMBOLS)
      USE MODULE_PREC(RAL_NLLS_ROOTS), ONLY: ROOTS_cubic

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: DTRS_initialize, DTRS_solve, DTRS_solve_main

!--------------------
!   P r e c i s i o n
!--------------------

!     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: history_max = 100
      INTEGER, PARAMETER :: max_degree = 3
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: point4 = 0.4_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: six = 6.0_wp
      REAL ( KIND = wp ), PARAMETER :: sixth = one / six
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: twentyfour = 24.0_wp
      REAL ( KIND = wp ), PARAMETER :: largest = HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: lower_default = - half * largest
      REAL ( KIND = wp ), PARAMETER :: upper_default = largest
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: teneps = ten * epsmch
      REAL ( KIND = wp ), PARAMETER :: roots_tol = teneps
      LOGICAL :: roots_debug = .FALSE.

!--------------------------
!  Derived type definitions
!--------------------------

!  - - - - - - - - - - - - - - - - - - - - - - -
!   control derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_control_type

!  unit for error messages

        INTEGER :: error = 6

!  unit for monitor output

        INTEGER :: out = 6

!  unit to write problem data into file problem_file

        INTEGER :: problem = 0

!  controls level of diagnostic output

        INTEGER :: print_level = 0

!  maximum degree of Taylor approximant allowed

        INTEGER :: taylor_max_degree = 3

!  any entry of H that is smaller than h_min * MAXVAL( H ) we be treated as zero

        REAL ( KIND = wp ) :: h_min = epsmch

!  any entry of C that is smaller than c_min * MAXVAL( C ) we be treated as zero

        REAL ( KIND = wp ) :: c_min = epsmch

!  lower and upper bounds on the multiplier, if known

        REAL ( KIND = wp ) :: lower = lower_default
        REAL ( KIND = wp ) :: upper = upper_default

!  stop when | ||x|| - radius | <=
!     max( stop_normal * radius, stop_absolute_normal )

        REAL ( KIND = wp ) :: stop_normal = epsmch
        REAL ( KIND = wp ) :: stop_absolute_normal = epsmch

!  is the solution is REQUIRED to lie on the boundary (i.e., is the constraint
!  an equality)?

        LOGICAL :: equality_problem= .FALSE.

!  name of file into which to write problem data

        CHARACTER ( LEN = 30 ) :: problem_file =                               &
         'trs_problem.data' // REPEAT( ' ', 14 )

!  all output lines will be prefixed by
!    prefix(2:LEN(TRIM(%prefix))-1)
!  where prefix contains the required string enclosed in quotes,
!  e.g. "string" or 'string'

        CHARACTER ( LEN = 30 ) :: prefix  = '""                            '

      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - - -
!   history derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_history_type

!  value of lambda

        REAL ( KIND = wp ) :: lambda = zero

!  corresponding value of ||x(lambda)||_M

        REAL ( KIND = wp ) :: x_norm = zero
      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - -
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DTRS_inform_type

!   reported return status:
!      0 the solution has been found
!     -3 n and/or Delta is not positive
!    -16 ill-conditioning has prevented furthr progress

        INTEGER :: status = 0

!  the number of (||x||_M,lambda) pairs in the history

        INTEGER :: len_history = 0

!  the value of the quadratic function

        REAL ( KIND = wp ) :: obj = HUGE( one )

!  the M-norm of x, ||x||_M

        REAL ( KIND = wp ) :: x_norm = zero

!  the Lagrange multiplier corresponding to the trust-region constraint

        REAL ( KIND = wp ) :: multiplier = zero

!  a lower bound max(0,-lambda_1), where lambda_1 is the left-most
!  eigenvalue of (H,M)

        REAL ( KIND = wp ) :: pole = zero

!  has the hard case occurred?

        LOGICAL :: hard_case = .FALSE.

!  history information

        TYPE ( DTRS_history_type ), DIMENSION( history_max ) :: history
      END TYPE

!  interface to BLAS: two norm

     INTERFACE NRM2

       FUNCTION SNRM2( n, X, incx )
       import :: lp
       REAL(KIND=lp) :: SNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL(KIND=lp), INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION SNRM2

       FUNCTION DNRM2( n, X, incx )
       import :: np
       REAL(KIND=np):: DNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL(KIND=np), INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION DNRM2

     END INTERFACE

    CONTAINS

!-*-*-*-*-*-*-  D T R S _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE DTRS_initialize( control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  .  Set initial values for the TRS control parameters  .
!
!  Arguments:
!  =========
!
!   control  a structure containing control information. See DTRS_control_type
!   inform   a structure containing information. See DRQS_inform_type
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------

      TYPE ( DTRS_CONTROL_TYPE ), INTENT( OUT ) :: control
      TYPE ( DTRS_inform_type ), INTENT( OUT ) :: inform

      inform%status = RAL_NLLS_ok

!  Set initial control parameter values

      control%stop_normal = epsmch ** 0.75
      control%stop_absolute_normal = epsmch ** 0.75

      RETURN

!  End of subroutine DTRS_initialize

      END SUBROUTINE DTRS_initialize

!-*-*-*-*-*-*-*-*-  D T R S _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*-

      SUBROUTINE DTRS_solve( n, radius, f, C, H, X, control, inform,           &
                             C_scale, H_scale )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the trust-region subproblem
!
!      minimize     q(x) = 1/2 <x, H x> + <c, x> + f
!      subject to   ||x||_2 <= radius  or ||x||_2 = radius
!
!  where H is diagonal, using a secular iteration
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   radius - the trust-region radius
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a vector of values for the diagonal matrix H
!
!   X - the required solution vector x
!
!   control - a structure containing control information. See DTRS_control_type
!
!   inform - a structure containing information. See DTRS_inform_type
!
!   C_scale, H_scale - workspace arrays
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: radius
      REAL ( KIND = wp ), INTENT( IN ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C, H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( DTRS_control_type ), INTENT( IN ) :: control
      TYPE ( DTRS_inform_type ), INTENT( INOUT ) :: inform
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: C_scale, H_scale

!  local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: scale_h, scale_c, f_scale, radius_scale
      TYPE ( DTRS_control_type ) :: control_scale

!  scale the problem to solve instead

!      minimize    q_s(x_s) = 1/2 <x_s, H_s x_s> + <c_s, x_s> + f_s
!      subject to    ||x_s||_2 <= radius_s  or ||x_s||_2 = radius_s

!  where H_s = H / s_h and c_s = c / s_c for scale factors s_h and s_c

!  This corresponds to
!    radius_s = ( s_h / s_c ) radius,
!    f_s = ( s_h / s_c^2 ) f
!  and the solution may be recovered as
!    x = ( s_c / s_h ) x_s
!    lambda = s_h lambda_s
!    q(x) = ( s_c^2 / s_ h ) q_s(x_s)

!write(6,"( A2, 5ES13.4E3 )" ) 'H', H
!write(6,"( A2, 5ES13.4E3 )" ) 'C', C

!  scale H by the largest H and remove relatively tiny H

     !scale_h = MAXVAL( ABS( H ) )
      scale_h = zero
      DO i = 1, n
        scale_h = MAX( scale_h, ABS( H( i ) ) )
      END DO

      IF ( scale_h > zero ) THEN
        DO i = 1, n
          IF ( ABS( H( i ) ) >= control%h_min * scale_h ) THEN
            H_scale( i ) = H( i ) / scale_h
          ELSE
            H_scale( i ) = zero
          END IF
        END DO
      ELSE
        scale_h = one
        H_scale = zero
      END IF

!  scale c by the largest c and remove relatively tiny c

     !scale_c = MAXVAL( ABS( C ) )
      scale_c = zero
      DO i = 1, n
        scale_c = MAX( scale_c, ABS( C( i ) ) )
      END DO

      IF ( scale_c > zero ) THEN
        DO i = 1, n
          IF ( ABS( C( i ) ) >= control%h_min * scale_c ) THEN
            C_scale( i ) = C( i ) / scale_c
          ELSE
            C_scale( i ) = zero
          END IF
        END DO
      ELSE
        scale_c = one
        C_scale = zero
      END IF

      radius_scale = ( scale_h / scale_c ) * radius
      f_scale = ( scale_h / scale_c ** 2 ) * f

      control_scale = control
      IF ( control_scale%lower /= lower_default )                              &
        control_scale%lower = control_scale%lower / scale_h
      IF ( control_scale%upper /= upper_default )                              &
        control_scale%upper = control_scale%upper / scale_h

!  solve the scaled problem

      CALL DTRS_solve_main( n, radius_scale, f_scale, C_scale, H_scale, X,     &
                            control_scale, inform )

!  unscale the solution, function value, multiplier and related values

      X = ( scale_c / scale_h ) * X
      inform%obj = ( scale_c ** 2 / scale_h ) * inform%obj
      inform%multiplier = scale_h * inform%multiplier
      inform%pole = scale_h * inform%pole
      DO i = 1, inform%len_history
        inform%history( i )%lambda = scale_h * inform%history( i )%lambda
        inform%history( i )%x_norm                                             &
          = ( scale_c / scale_h ) * inform%history( i )%x_norm
      END DO

!  End of subroutine DTRS_solve

      END SUBROUTINE DTRS_solve

!-*-*-*-*-*-*-*-*-  D T R S _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*-

      SUBROUTINE DTRS_solve_main( n, radius, f, C, H, X, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the trust-region subproblem
!
!      minimize     1/2 <x, H x> + <c, x> + f
!      subject to    ||x||_2 <= radius  or ||x||_2 = radius
!
!  where H is diagonal, using a secular iteration
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   radius - the trust-region radius
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a vector of values for the diagonal matrix H
!
!   X - the required solution vector x
!
!   control - a structure containing control information. See DTRS_control_type
!
!   inform - a structure containing information. See DTRS_inform_type
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: radius
      REAL ( KIND = wp ), INTENT( IN ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C, H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( DTRS_control_type ), INTENT( IN ) :: control
      TYPE ( DTRS_inform_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, it, out, nroots, print_level
      INTEGER :: max_order, n_lambda, i_hard, n_sing
      REAL ( KIND = wp ) :: lambda, lambda_l, lambda_u, delta_lambda
      REAL ( KIND = wp ) :: alpha, utx, distx, x_big
      REAL ( KIND = wp ) :: c_norm, c_norm_over_radius, v_norm2, w_norm2, h_norm
      REAL ( KIND = wp ) :: beta, z_norm2, root1, root2, root3
      REAL ( KIND = wp ) :: lambda_min, lambda_max, lambda_plus, c2
      REAL ( KIND = wp ) :: a_0, a_1, a_2, a_3, a_max
      REAL ( KIND = wp ), DIMENSION( 3 ) :: lambda_new
      INTEGER, DIMENSION( 1 ) :: mloc
      REAL ( KIND = wp ), DIMENSION( 0 : max_degree ) :: x_norm2, pi_beta
      LOGICAL :: printi, printt, printd, problem_file_exists
      CHARACTER ( LEN = 1 ) :: region

!  prefix for all output

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      IF ( LEN( TRIM( control%prefix ) ) > 2 )                                 &
        prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!write(6,"( A2, 5ES13.4E3 )" ) 'H', H
!write(6,"( A2, 5ES13.4E3 )" ) 'C', C
!write(6,"( A, ES13.4E3 )" ) 'radius', radius

!  output problem data

      IF ( control%problem > 0 ) THEN
        INQUIRE( FILE = control%problem_file, EXIST = problem_file_exists )
        IF ( problem_file_exists ) THEN
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'OLD' )
          REWIND control%problem
        ELSE
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'NEW' )
        END IF
        WRITE( control%problem, * ) n, COUNT( C( : n ) /= zero ),              &
          COUNT( H( : n ) /= zero )
        WRITE( control%problem, * ) radius, f
        DO i = 1, n
          IF ( C( i ) /= zero ) WRITE( control%problem, * ) i, C( i )
        END DO
        DO i = 1, n
          IF ( H( i ) /= zero ) WRITE( control%problem, * ) i, i, H( i )
        END DO
        CLOSE( control%problem )
      END IF

!  set initial values

      X = zero ; inform%x_norm = zero ; inform%obj = f
      inform%hard_case = .FALSE.
      delta_lambda = zero

!  record desired output level

      out = control%out
      print_level = control%print_level
      printi = out > 0 .AND. print_level > 0
      printt = out > 0 .AND. print_level > 1
      printd = out > 0 .AND. print_level > 2

!  check for n < 0 or radius < 0

      IF ( n < 0 .or. radius < 0 ) THEN
        inform%status = RAL_NLLS_error_restrictions
        RETURN
      END IF

!  compute the two-norm of c and the extreme eigenvalues of H

      c_norm = TWO_NORM( C )
      lambda_min = MINVAL( H( : n ) )
      lambda_max = MAXVAL( H( : n ) )
      h_norm = zero
      DO i = 1, n
        h_norm = MAX( h_norm, ABS( H( i ) ) )
      END DO
      IF ( printt ) WRITE( out, "( A, ' ||c|| = ', ES10.4, ', ||H|| = ',       &
     &                             ES10.4, ', lambda_min = ', ES11.4 )" )      &
          prefix, c_norm, h_norm, lambda_min

      region = 'L'
      IF ( printt )                                                            &
        WRITE( out, "( A, 4X, 28( '-' ), ' phase two ', 28( '-' ) )" ) prefix
      IF ( printi ) WRITE( out, 2030 ) prefix

!  check for the trivial case

      IF ( c_norm == zero .AND. lambda_min >= zero ) THEN
        IF (  control%equality_problem ) THEN
          DO i = 1, n
            IF ( H( i ) == lambda_min ) THEN
              i_hard = i
              EXIT
            END IF
          END DO
          X( i_hard ) = one / radius
          inform%x_norm = radius
          inform%obj = f + lambda_min * radius ** 2
          lambda = - lambda_min
        ELSE
          lambda = zero
        END IF
        IF ( printi ) THEN
          WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,             &
          0, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
          WRITE( out, "( A,                                                    &
       &    ' Normal stopping criteria satisfied' )" ) prefix
        END IF
        inform%status = RAL_NLLS_ok
        GO TO 900
      END IF

!  construct values lambda_l and lambda_u for which lambda_l <= lambda_optimal
!   <= lambda_u, and ensure that all iterates satisfy lambda_l <= lambda
!   <= lambda_u

      c_norm_over_radius = c_norm / radius
      IF ( control%equality_problem ) THEN
        lambda_l = MAX( control%lower,  - lambda_min,                          &
                        c_norm_over_radius - lambda_max )
        lambda_u = MIN( control%upper,                                         &
                        c_norm_over_radius - lambda_min )
      ELSE
        lambda_l = MAX( control%lower, zero, - lambda_min,                     &
                        c_norm_over_radius - lambda_max )
        lambda_u = MIN( control%upper,                                         &
                        MAX( zero, c_norm_over_radius - lambda_min ) )
      END IF
      lambda = lambda_l

!  check for the "hard case"

      IF ( lambda == - lambda_min ) THEN
        c2 = zero
        inform%hard_case = .TRUE.
        DO i = 1, n
          IF ( H( i ) == lambda_min ) THEN
!           IF ( ABS( C( i ) ) > epsmch ) THEN
            IF ( ABS( C( i ) ) > epsmch * c_norm ) THEN
              inform%hard_case = .FALSE.
              c2 = c2 + C( i ) ** 2
            ELSE
              i_hard = i
            END IF
          END IF
        END DO

!  the hard case may occur

        IF ( inform%hard_case ) THEN
          DO i = 1, n
            IF ( H( i ) /= lambda_min ) THEN
              X( i )  = - C( i ) / ( H( i ) + lambda )
            ELSE
              X( i ) = zero
            END IF
          END DO
          inform%x_norm = TWO_NORM( X )

!  the hard case does occur

          IF ( inform%x_norm <= radius ) THEN
            IF ( inform%x_norm < radius ) THEN

!  compute the step alpha so that X + alpha E_i_hard lies on the trust-region
!  boundary and gives the smaller value of q

              utx = X( i_hard ) / radius
              distx = ( radius - inform%x_norm ) *                             &
                ( ( radius + inform%x_norm ) / radius )
              alpha = sign( distx / ( abs( utx ) +                             &
                            sqrt( utx ** 2 + distx / radius ) ), utx )

!  record the optimal values

              X( i_hard ) = X( i_hard ) + alpha
            END IF
            inform%x_norm = TWO_NORM( X )
            inform%obj =                                                       &
                f + half * ( DOT_PRODUCT( C, X ) - lambda * radius ** 2 )
            IF ( printi ) THEN
              WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,         &
              0, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
              WRITE( out, "( A,                                                &
           &    ' Hard-case stopping criteria satisfied' )" ) prefix
            END IF
            inform%status = RAL_NLLS_ok
            GO TO 900

!  the hard case didn't occur after all

          ELSE
            inform%hard_case = .FALSE.

!  compute the first derivative of ||x|(lambda)||^2 - radius^2

            w_norm2 = zero
            DO i = 1, n
              IF ( H( i ) /= lambda_min )                                      &
                w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
            END DO
            x_norm2( 1 ) = - two * w_norm2

!  compute the Newton correction

            lambda =                                                           &
              lambda - ( inform%x_norm ** 2 - radius ** 2 ) / x_norm2( 1 )
            lambda_l = MAX( lambda_l, lambda )
          END IF

!  there is a singularity at lambda. Compute the point for which the
!  sum of squares of the singular terms is equal to radius^2

        ELSE
!         lambda = lambda + SQRT( c2 ) / radius
          lambda = lambda + MAX( SQRT( c2 ) / radius, lambda * epsmch )
          lambda_l = MAX( lambda_l, lambda )
        END IF
      END IF

!  the iterates will all be in the L region. Prepare for the main loop

      it = 0
      max_order = MAX( 1, MIN( max_degree, control%taylor_max_degree ) )

!  start the main loop

      DO
        it = it + 1
!       IF ( it > 10 ) stop

!  if H(lambda) is positive definite, solve  H(lambda) x = - c

        n_sing = COUNT( H( : n ) + lambda <= zero )
        IF ( n_sing == 0 ) THEN
          DO i = 1, n
            X( i )  = - C( i ) / ( H( i ) + lambda )
          END DO
        ELSE
          x_big = radius / SQRT( REAL( n_sing, wp ) )
          DO i = 1, n
            IF ( H( i ) + lambda > zero ) THEN
              X( i )  = - C( i ) / ( H( i ) + lambda )
            ELSE
              X( i )  = SIGN( x_big, - C( i ) )
            END IF
          END DO
        END IF
!write(6,"( A2, 5ES13.4E3 )" ) 'X', X

!  compute the two-norm of x

        inform%x_norm = TWO_NORM( X )
        x_norm2( 0 ) = inform%x_norm ** 2

!  if the Newton step lies within the trust region, exit

        IF ( lambda == zero .AND. inform%x_norm <= radius ) THEN
          inform%obj = f + half * DOT_PRODUCT( C, X )
          inform%status = RAL_NLLS_ok
          region = 'L'
          IF ( printi ) THEN
            WRITE( out, "( A, A2, I4, 2ES22.15 )" ) prefix,                    &
              region, it, inform%x_norm - radius, lambda
            WRITE( out, "( A, ' Interior stopping criteria satisfied')" ) prefix
          END IF
          GO TO 900
        END IF

!  the current estimate gives a good approximation to the required root

        IF ( ABS( inform%x_norm - radius ) <=                                  &
               MAX( control%stop_normal * radius,                              &
                    control%stop_absolute_normal ) ) THEN
          IF ( inform%x_norm > radius ) THEN
            lambda_l = MAX( lambda_l, lambda )
          ELSE
            region = 'G'
            lambda_u = MIN( lambda_u, lambda )
          END IF
          IF ( printt .AND. it > 1 ) WRITE( out, 2030 ) prefix
          IF ( printi ) THEN
            WRITE( out, "( A, A2, I4, 3ES22.15 )" )  prefix, region,           &
            it, ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
            WRITE( out, "( A,                                                  &
         &    ' Normal stopping criteria satisfied' )" ) prefix
          END IF
          inform%status = RAL_NLLS_ok
          EXIT
        END IF

        lambda_l = MAX( lambda_l, lambda )

!  debug printing

        IF ( printd ) THEN
          WRITE( out, "( A, 8X, 'lambda', 13X, 'x_norm', 15X, 'radius' )" )    &
            prefix
          WRITE( out, "( A, 3ES20.12 )") prefix, lambda, inform%x_norm, radius
        ELSE IF ( printi ) THEN
          IF ( printt .AND. it > 1 ) WRITE( out, 2030 ) prefix
          WRITE( out, "( A, A2, I4, 3ES22.15 )" ) prefix, region, it,          &
            ABS( inform%x_norm - radius ), lambda, ABS( delta_lambda )
        END IF

!  record, for the future, values of lambda which give small ||x||

        IF ( inform%len_history < history_max ) THEN
          inform%len_history = inform%len_history + 1
          inform%history( inform%len_history )%lambda = lambda
          inform%history( inform%len_history )%x_norm = inform%x_norm
        END IF

!  a lambda in L has been found. It is now simply a matter of applying
!  a variety of Taylor-series-based methods starting from this lambda

!  precaution against rounding producing lambda outside L

        IF ( lambda > lambda_u ) THEN
          inform%status = RAL_NLLS_error_ill_conditioned
          EXIT
        END IF

!  compute first derivatives of x^T M x

!  form ||w||^2 = x^T H^-1(lambda) x

        w_norm2 = zero
        DO i = 1, n
          w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
        END DO

!  compute the first derivative of x_norm2 = x^T M x

        x_norm2( 1 ) = - two * w_norm2

!  compute pi_beta = ||x||^beta and its first derivative when beta = - 1

        beta = - one
        CALL DTRS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )

!  compute the Newton correction (for beta = - 1)

        delta_lambda = - ( pi_beta( 0 ) - ( radius ) ** beta ) / pi_beta( 1 )

        n_lambda = 1
        lambda_new( n_lambda ) = lambda + delta_lambda

        IF ( max_order >= 3 ) THEN

!  compute the second derivative of x^T x

          z_norm2 = zero
          DO i = 1, n
            z_norm2 = z_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 4
          END DO
          x_norm2( 2 ) = six * z_norm2

!  compute the third derivatives of x^T x

          v_norm2 = zero
          DO i = 1, n
            v_norm2 = v_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 5
          END DO
          x_norm2( 3 ) = - twentyfour * v_norm2

!  compute pi_beta = ||x||^beta and its derivatives when beta = 2

          beta = two
          CALL DTRS_pi_derivs( max_order, beta, x_norm2( : max_order ),        &
                               pi_beta( : max_order ) )

!  compute the "cubic Taylor approximaton" step (beta = 2)

          a_0 = pi_beta( 0 ) - ( radius ) ** beta
          a_1 = pi_beta( 1 )
          a_2 = half * pi_beta( 2 )
          a_3 = sixth * pi_beta( 3 )
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max
            a_2 = a_2 / a_max ; a_3 = a_3 / a_max
          END IF
          CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,             &
                            root1, root2, root3, roots_debug )
          n_lambda = n_lambda + 1
          IF ( nroots == 3 ) THEN
            lambda_new( n_lambda ) = lambda + root3
          ELSE
            lambda_new( n_lambda ) = lambda + root1
          END IF

!  compute pi_beta = ||x||^beta and its derivatives when beta = - 0.4

          beta = - point4
          CALL DTRS_pi_derivs( max_order, beta, x_norm2( : max_order ),        &
                               pi_beta( : max_order ) )

!  compute the "cubic Taylor approximaton" step (beta = - 0.4)

          a_0 = pi_beta( 0 ) - ( radius ) ** beta
          a_1 = pi_beta( 1 )
          a_2 = half * pi_beta( 2 )
          a_3 = sixth * pi_beta( 3 )
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max
            a_2 = a_2 / a_max ; a_3 = a_3 / a_max
          END IF
          CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,             &
                            root1, root2, root3, roots_debug )
          n_lambda = n_lambda + 1
          IF ( nroots == 3 ) THEN
            lambda_new( n_lambda ) = lambda + root3
          ELSE
            lambda_new( n_lambda ) = lambda + root1
          END IF
        END IF

!  record all of the estimates of the optimal lambda

        IF ( printd ) THEN
          MLOC(:) = MAXLOC( lambda_new( : n_lambda ) )
          WRITE( out, "( A, ' lambda_t (', I1, ')', 3ES20.13 )" )              &
            prefix, MLOC(1), lambda_new( : MIN( 3, n_lambda ) )
        END IF

!  compute the best Taylor improvement

        lambda_plus = MAXVAL( lambda_new( : n_lambda ) )
        delta_lambda = lambda_plus - lambda
        lambda = lambda_plus

!  improve the lower bound if possible

        lambda_l = MAX( lambda_l, lambda_plus )

!  check that the best Taylor improvement is significant

        IF ( ABS( delta_lambda ) < epsmch * MAX( one, ABS( lambda ) ) ) THEN
          IF ( printi ) WRITE( out, "( A, ' normal exit with no ',             &
         &                     'significant Taylor improvement' )" ) prefix
          inform%status = RAL_NLLS_ok
          EXIT
        END IF

!  End of main iteration loop

      END DO

!  Record the optimal obective value

      inform%obj = f + half * ( DOT_PRODUCT( C, X ) - lambda * x_norm2( 0 ) )

!  ----
!  Exit
!  ----

 900  CONTINUE
      inform%multiplier = lambda
      inform%pole = MAX( zero, - lambda_min )
      RETURN

! Non-executable statements

 2030 FORMAT( A, '    it     ||x||-radius             lambda ',                &
                 '              d_lambda' )

!  End of subroutine DTRS_solve_main

      END SUBROUTINE DTRS_solve_main

!-*-*-*-*-*-*-  D T R S _ P I _ D E R I V S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE DTRS_pi_derivs( max_order, beta, x_norm2, pi_beta )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute pi_beta = ||x||^beta and its derivatives
!
!  Arguments:
!  =========
!
!  Input -
!   max_order - maximum order of derivative
!   beta - power
!   x_norm2 - (0) value of ||x||^2,
!             (i) ith derivative of ||x||^2, i = 1, max_order
!  Output -
!   pi_beta - (0) value of ||x||^beta,
!             (i) ith derivative of ||x||^beta, i = 1, max_order
!
!  Extracted wholesale from module RAL_NLLS_RQS
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: max_order
      REAL ( KIND = wp ), INTENT( IN ) :: beta, x_norm2( 0 : max_order )
      REAL ( KIND = wp ), INTENT( OUT ) :: pi_beta( 0 : max_order )

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      REAL ( KIND = wp ) :: hbeta

      hbeta = half * beta
      pi_beta( 0 ) = x_norm2( 0 ) ** hbeta
      pi_beta( 1 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - one ) ) * x_norm2( 1 )
      IF ( max_order == 1 ) RETURN
      pi_beta( 2 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - two ) ) *             &
        ( ( hbeta - one ) * x_norm2( 1 ) ** 2 + x_norm2( 0 ) * x_norm2( 2 ) )
      IF ( max_order == 2 ) RETURN
      pi_beta( 3 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - three ) ) *           &
        ( x_norm2( 3 ) * x_norm2( 0 ) ** 2 + ( hbeta - one ) *                 &
          ( three * x_norm2( 0 ) * x_norm2( 1 ) * x_norm2( 2 ) +               &
            ( hbeta - two ) * x_norm2( 1 ) ** 3 ) )

      RETURN

!  End of subroutine DTRS_pi_derivs

      END SUBROUTINE DTRS_pi_derivs

!-*-*-*-*-  G A L A H A D   T W O  _ N O R M   F U N C T I O N   -*-*-*-*-

       FUNCTION TWO_NORM( X )

!  Compute the l_2 norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: TWO_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         TWO_NORM = NRM2( n, X, 1 )
       ELSE
         TWO_NORM = zero
       END IF
       RETURN

!  End of function TWO_NORM

       END FUNCTION TWO_NORM

!-*-*-*-*-*-  End of R A L _ N L L S _ D T R S  double  M O D U L E  *-*-*-*-*-

   END MODULE MODULE_PREC(RAL_NLLS_DTRS)



! THIS VERSION: RAL_NLLS 1.2 - 04/03/2020 AT 13:15 GMT.

!-*-*-*-*-*-*-*-  R A L _ N L L S _ R Q S  double  M O D U L E  *-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package RQS, March 4th, 2016

!  For full documentation, see
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE MODULE_PREC(RAL_NLLS_DRQS)

!       --------------------------------------------
!      |                                            |
!      | Solve the reguarized quadratic subproblem  |
!      |                                            |
!      |    minimize     1/2 <x, H x> + <c, x> + f  |
!      |                   + (sigma/p) ||x||_2^p    |
!      |                                            |
!      | where H is diagonal                        |
!      |                                            |
!       --------------------------------------------

      USE MODULE_PREC(ral_nlls_types), Only: lp, np, wp
      USE MODULE_PREC(RAL_NLLS_SYMBOLS)
      USE MODULE_PREC(RAL_NLLS_ROOTS)

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: DRQS_initialize, DRQS_solve, DRQS_solve_main

!--------------------
!   P r e c i s i o n
!--------------------

!     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: history_max = 100
      INTEGER, PARAMETER :: max_degree = 3
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: point4 = 0.4_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: six = 6.0_wp
      REAL ( KIND = wp ), PARAMETER :: sixth = one / six
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: twentyfour = 24.0_wp
      REAL ( KIND = wp ), PARAMETER :: largest = HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: lower_default = - half * largest
      REAL ( KIND = wp ), PARAMETER :: upper_default = largest
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: teneps = ten * epsmch
      REAL ( KIND = wp ), PARAMETER :: roots_tol = teneps
      LOGICAL :: roots_debug = .FALSE.

!--------------------------
!  Derived type definitions
!--------------------------

!  - - - - - - - - - - - - - - - - - - - - - - -
!   control derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DRQS_control_type

!  unit for error messages

        INTEGER :: error = 6

!  unit for monitor output

        INTEGER :: out = 6

!  unit to write problem data into file problem_file

        INTEGER :: problem = 0

!  controls level of diagnostic output

        INTEGER :: print_level = 0

!   at most maxit inner iterations are allowed

        INTEGER :: maxit = 1000

!  maximum degree of Taylor approximant allowed

        INTEGER :: taylor_max_degree = 3

!  any entry of H that is smaller than h_min * MAXVAL( H ) we be treated as zero

        REAL ( KIND = wp ) :: h_min = epsmch ** 2

!  any entry of C that is smaller than c_min * MAXVAL( C ) we be treated as zero

        REAL ( KIND = wp ) :: c_min = epsmch ** 2

!  lower and upper bounds on the multiplier, if known

        REAL ( KIND = wp ) :: lower = - half * HUGE( one )
        REAL ( KIND = wp ) :: upper =  HUGE( one )

!  stop when | ||x|| - (multiplier/sigma)^(1/(p-2)) | <=
!    max( stop_normal * max( ||x||, (multiplier/sigma)^(1/(p-2)) ),
!         stop_absolute_normal )

        REAL ( KIND = wp ) :: stop_normal = epsmch
        REAL ( KIND = wp ) :: stop_absolute_normal = epsmch

!  name of file into which to write problem data

        CHARACTER ( LEN = 30 ) :: problem_file =                               &
         'rqs_problem.data' // REPEAT( ' ', 14 )

!  all output lines will be prefixed by
!    prefix(2:LEN(TRIM(%prefix))-1)
!  where prefix contains the required string enclosed in quotes,
!  e.g. "string" or 'string'

        CHARACTER ( LEN = 30 ) :: prefix  = '""                            '
      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - - -
!   history derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DRQS_history_type

!  value of lambda

        REAL ( KIND = wp ) :: lambda = zero

!  corresponding value of ||x(lambda)||_M

        REAL ( KIND = wp ) :: x_norm = zero
      END TYPE

!  - - - - - - - - - - - - - - - - - - - - - - -
!   inform derived type with component defaults
!  - - - - - - - - - - - - - - - - - - - - - - -

      TYPE, PUBLIC :: DRQS_inform_type

!   reported return status:
!      0 the solution has been found
!     -3 n and/or sigma is not positive and/or p < 2
!     -7 the problem is unbounded from below when p = 2
!    -16 ill-conditioning has prevented furthr progress

        INTEGER :: status = 0

!  STAT value after allocate failure

        INTEGER :: alloc_status = 0

!  the number of (||x||_M,lambda) pairs in the history

        INTEGER :: len_history = 0

!  the value of the quadratic function

        REAL ( KIND = wp ) :: obj = HUGE( one )

!  the value of the regularized quadratic function

        REAL ( KIND = wp ) :: obj_regularized = HUGE( one )

!  the M-norm of x, ||x||_M

        REAL ( KIND = wp ) :: x_norm = zero

!  the Lagrange multiplier corresponding to the regularization

        REAL ( KIND = wp ) :: multiplier = zero

!  a lower bound max(0,-lambda_1), where lambda_1 is the left-most
!  eigenvalue of (H,M)

        REAL ( KIND = wp ) :: pole = zero

!  has the hard case occurred?

        LOGICAL :: hard_case = .FALSE.

!  history information

        TYPE ( DRQS_history_type ), DIMENSION( history_max ) :: history
      END TYPE

!  interface to BLAS: two norm

     INTERFACE NRM2

       FUNCTION SNRM2( n, X, incx )
       import :: lp
       REAL(KIND=lp) :: SNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL(KIND=lp), INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION SNRM2

       FUNCTION DNRM2( n, X, incx )
       import :: np
       REAL(KIND=np):: DNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL(KIND=np), INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION DNRM2

     END INTERFACE

    CONTAINS

!-*-*-*-*-*-*-  D R Q S _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE DRQS_initialize( control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  .  Set initial values for the DRQS control parameters  .
!
!  Arguments:
!  =========
!
!   control  a structure containing control information. See DRQS_control_type
!   inform   a structure containing information. See DRQS_inform_type
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!---------------------------------
!   D u m m y   A r g u m e n t s
!---------------------------------

      TYPE ( DRQS_CONTROL_TYPE ), INTENT( OUT ) :: control
      TYPE ( DRQS_inform_type ), INTENT( OUT ) :: inform

      inform%status = RAL_NLLS_ok

!  Set initial control parameter values

      control%stop_normal = epsmch ** 0.75
      control%stop_absolute_normal = epsmch ** 0.75

      RETURN

!  End of subroutine DRQS_initialize

      END SUBROUTINE DRQS_initialize

!-*-*-*-*-*-*-*-*-  D R Q S _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*-*-

      SUBROUTINE DRQS_solve( n, p, sigma, f, C, H, X, control, inform,         &
                             C_scale, H_scale )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the reguarized quadratic subproblem
!
!      minimize     1/2 <x, H x> + <c, x> + f + (sigma/p) ||x||_2^p
!
!  where H is diagonal using a secular iteration
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   p - the order of the regularization
!
!   sigma - the regularization weight
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a vector of values for the diagonal matrix H
!
!   X - the required solution vector x
!
!   control - a structure containing control information. See DRQS_control_type
!
!   inform - a structure containing information. See DRQS_inform_type
!
!   C_scale, H_scale - workspace arrays
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: p, sigma, f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C, H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( DRQS_control_type ), INTENT( IN ) :: control
      TYPE ( DRQS_inform_type ), INTENT( INOUT ) :: inform
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: C_scale, H_scale

!  local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: scale_h, scale_c, f_scale, sigma_scale
      TYPE ( DRQS_control_type ) :: control_scale

!  scale the problem to solve instead

!  minimize    q_s(x_s) = 1/2 <x_s, H_s x_s> + <c_s, x_s> + f_s
!                + (sigma_s/p) ||x_s||_2^p

!  where H_s = H / s_h and c_s = c / s_c for scale factors s_h and s_c

!  This corresponds to
!    sigma_s = ( s_h / s_c^2 ) ( s_c / s_h )^p,
!    f_s = ( s_h / s_c^2 ) f
!  and the solution may be recovered as
!    x = ( s_c / s_h ) x_s
!    lambda = s_h lambda_s
!    q(x) = ( s_c^2 / s_ h ) q_s(x_s)

!write(6,"( ' sigma ', ES13.4E3 )" ) sigma
!write(6,"( A2, 5ES13.4E3 )" ) 'H', H
!write(6,"( A2, 5ES13.4E3 )" ) 'C', C

!  scale H by the largest H and remove relatively tiny H

     !scale_h = MAXVAL( ABS( H ) )
      scale_h = zero
      DO i = 1, n
        scale_h = MAX( scale_h, ABS( H( i ) ) )
      END DO
      IF ( scale_h > zero ) THEN
        DO i = 1, n
          IF ( ABS( H( i ) ) >= control%h_min * scale_h ) THEN
            H_scale( i ) = H( i ) / scale_h
          ELSE
            H_scale( i ) = zero
          END IF
        END DO
      ELSE
        scale_h = one
        H_scale = zero
      END IF

!  scale c by the largest c and remove relatively tiny c

     !scale_c = MAXVAL( ABS( C ) )
      scale_c = zero
      DO i = 1, n
        scale_c = MAX( scale_c, ABS( C( i ) ) )
      END DO
      IF ( scale_c > zero ) THEN
        DO i = 1, n
          IF ( ABS( C( i ) ) >= control%h_min * scale_c ) THEN
            C_scale( i ) = C( i ) / scale_c
          ELSE
            C_scale( i ) = zero
          END IF
        END DO
      ELSE
        scale_c = one
        C_scale = zero
      END IF

      sigma_scale                                                              &
        = sigma * ( scale_h / scale_c ** 2 ) * ( scale_c / scale_h ) ** p
      f_scale = f * ( scale_h / scale_c ** 2 )

      control_scale = control
      IF ( control_scale%lower /= lower_default )                              &
        control_scale%lower = control_scale%lower / scale_h
      IF ( control_scale%upper /= upper_default )                              &
        control_scale%upper = control_scale%upper / scale_h

      CALL DRQS_solve_main( n, p, sigma_scale, f_scale, C_scale, H_scale, X,   &
                            control_scale, inform )

!  unscale the solution, function value, multiplier and related values

      X = ( scale_c / scale_h ) * X
      inform%obj = ( scale_c ** 2 / scale_h ) * inform%obj
      inform%obj_regularized                                                   &
        = ( scale_c ** 2 / scale_h ) * inform%obj_regularized
      inform%multiplier = scale_h * inform%multiplier
      inform%pole = scale_h * inform%pole
      DO i = 1, inform%len_history
        inform%history( i )%lambda = scale_h * inform%history( i )%lambda
        inform%history( i )%x_norm                                             &
          = ( scale_c / scale_h ) * inform%history( i )%x_norm
      END DO

      RETURN

!  End of subroutine DRQS_solve

      END SUBROUTINE DRQS_solve

!-*-*-*-*-*-*-  D R Q S _ S O L V E _ M A I N   S U B R O U T I N E  -*-*-*-*-*-

      SUBROUTINE DRQS_solve_main( n, p, sigma, f, C, H, X, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the reguarized quadratic subproblem
!
!      minimize     1/2 <x, H x> + <c, x> + f + (sigma/p) ||x||_2^p
!
!  where H is diagonal using a secular iteration
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   p - the order of the regularization
!
!   sigma - the regularization weight
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a vector of values for the diagonal matrix H
!
!   X - the required solution vector x
!
!   control - a structure containing control information. See DRQS_control_type
!
!   inform - a structure containing information. See DRQS_inform_type
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: p, sigma, f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C, H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( DRQS_control_type ), INTENT( IN ) :: control
      TYPE ( DRQS_inform_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, it, out, nroots, print_level, max_order, n_lambda, i_hard
      REAL ( KIND = wp ) :: lambda, lambda_l, lambda_u, delta_lambda, target
      REAL ( KIND = wp ) :: alpha, utx, distx, c_norm, v_norm2, w_norm2, h_norm
      REAL ( KIND = wp ) :: beta, z_norm2, pm2, oopm2, oos, oos2
      REAL ( KIND = wp ) :: lambda_min, lambda_max, lambda_plus, topm2
      REAL ( KIND = wp ) :: a_0, a_1, a_2, a_3, a_max, c2
      REAL ( KIND = wp ), DIMENSION( 4 ) :: roots
      REAL ( KIND = wp ), DIMENSION( 9 ) :: lambda_new
      INTEGER, DIMENSION( 1 ) :: mloc
      REAL ( KIND = wp ), DIMENSION( 0 : max_degree ) :: x_norm2
      REAL ( KIND = wp ), DIMENSION( 0 : max_degree ) :: pi_beta, theta_beta
      LOGICAL :: printi, printt, printd, problem_file_exists
      CHARACTER ( LEN = 1 ) :: region

      INTEGER :: ii( 1 ), j
      REAL ( KIND = wp ) :: a_4

!  prefix for all output

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      IF ( LEN( TRIM( control%prefix ) ) > 2 )                                 &
        prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  output problem data

       IF ( control%problem > 0 ) THEN
        INQUIRE( FILE = control%problem_file, EXIST = problem_file_exists )
        IF ( problem_file_exists ) THEN
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'OLD' )
          REWIND control%problem
        ELSE
          OPEN( control%problem, FILE = control%problem_file,                  &
                FORM = 'FORMATTED', STATUS = 'NEW' )
        END IF
        WRITE( control%problem, * ) n, COUNT( C( : n ) /= zero ),              &
          COUNT( H( : n ) /= zero )
        WRITE( control%problem, * ) p, sigma, f
        DO i = 1, n
          IF ( C( i ) /= zero ) WRITE( control%problem, * ) i, C( i )
        END DO
        DO i = 1, n
          IF ( H( i ) /= zero ) WRITE( control%problem, * ) i, i, H( i )
        END DO
        CLOSE( control%problem )
      END IF

!write(6,"( ' sigma ', ES13.4E3 )" ) sigma
!write(6,"( A2, 5ES13.4E3 )" ) 'H', H
!write(6,"( A2, 5ES13.4E3 )" ) 'C', C

!  set initial values

      X = zero ; inform%x_norm = zero
      inform%obj = f ; inform%obj_regularized = f

      inform%hard_case = .FALSE.
      delta_lambda = zero

!  record desired output level

      out = control%out
      print_level = control%print_level
      printi = out > 0 .AND. print_level > 0
      printt = out > 0 .AND. print_level > 1
      printd = out > 0 .AND. print_level > 2

!  check for n < 0 or sigma < 0

      IF ( n < 0 .or. sigma < 0 ) THEN
        inform%status = RAL_NLLS_error_restrictions
        RETURN
      END IF

!  compute the two-norm of c and the extreme eigenvalues of H

      c_norm = TWO_NORM( C )
      lambda_min = MINVAL( H( : n ) )
      lambda_max = MAXVAL( H( : n ) )
      h_norm = zero
      DO i = 1, n
        h_norm = MAX( h_norm, ABS( H( i ) ) )
      END DO

      IF ( printi ) WRITE( out, "( A, ' ||c|| = ', ES10.4, ', ||H|| = ',       &
     &                             ES10.4, ', lambda_min = ', ES11.4 )" )      &
          prefix, c_norm, h_norm, lambda_min

      region = 'L'
      IF ( printt )                                                            &
        WRITE( out, "( A, 4X, 28( '-' ), ' phase two ', 28( '-' ) )" ) prefix
      IF ( printi ) WRITE( out, 2030 ) prefix

!  check for the trivial cases: ||c|| = 0 & H positive semi-definite

      IF ( c_norm == zero .AND. lambda_min >= zero ) THEN
        lambda = zero ; target = zero
        IF ( printi ) THEN
          WRITE( out, "( A, A2, I4, 1X, 3ES22.15 )" ) prefix, region,          &
            it, inform%x_norm - target, lambda, ABS( delta_lambda )
          WRITE( out, "( A,                                                    &
      &    ' Normal stopping criteria satisfied' )" ) prefix
        END IF
        inform%status = RAL_NLLS_ok
        GO TO 900
      END IF

!  p = 2

      IF ( p == two ) THEN
        IF ( lambda_min + sigma > zero ) THEN
          X = - C / ( H + sigma )
        ELSE IF ( lambda_min + sigma < zero ) THEN
          inform%status = RAL_NLLS_error_unbounded
          lambda = zero
          GO TO 900
        ELSE
          DO i = 1, n
            IF ( H( i ) + sigma <= zero .AND. C( i ) /= zero ) THEN
              inform%status = RAL_NLLS_error_unbounded
              lambda = zero
              GO TO 900
            ELSE
              X( i ) = - C( i ) / ( H( i ) + sigma )
            END IF
          END DO
        END IF
        lambda = sigma ; target = sigma
        inform%x_norm = TWO_NORM( X )
        inform%obj_regularized = f + half * DOT_PRODUCT( C, X )
        inform%obj = inform%obj_regularized - half * sigma * inform%x_norm ** 2
        inform%status = RAL_NLLS_ok
        GO TO 900
      END IF

!  reccord useful constants

      oos = one / sigma ; oos2 = oos * oos
      pm2 = p - two ; oopm2 = one / pm2 ; topm2 = two / pm2

!  construct values lambda_l and lambda_u for which lambda_l <= lambda_optimal
!   <= lambda_u, and ensure that all iterates satisfy lambda_l <= lambda
!   <= lambda_u

      lambda_l =                                                               &
        MAX( control%lower, zero, - lambda_min,                                &
             DRQS_lambda_root( lambda_max, c_norm * sigma ** oopm2, oopm2 ) )
      lambda_u =                                                               &
        MIN( control%upper, MAX( zero,                                         &
             DRQS_lambda_root( lambda_min, c_norm * sigma ** oopm2, oopm2 ) ) )
      lambda = lambda_l
!     write( 6,*) ' initial lambda ', lambda

!  find a better starting point for the p = 3 case

      IF ( p == 3 ) THEN
        DO i = 1, n
          a_0 = - sigma * ABS( C( i ) )
          a_1 = H( i )
          a_2 = one
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max ; a_2 = a_2 / a_max
          END IF
          CALL ROOTS_quadratic( a_0, a_1, a_2, roots_tol, nroots,              &
                                roots( 1 ), roots( 2 ), roots_debug )
          lambda = MAX( lambda, roots( 2 ) )
        END DO
!       write( 6,*) ' improved lambda ', lambda
      END IF

!write(6,*) ' lambda_l, lambda_u ', lambda_l, lambda_u

!lambda = SQRT( ABS( C( 1 ) * sigma ) )
!lambda = 1.0D-22

!  check for the "hard case"

      IF ( lambda == - lambda_min ) THEN
        c2 = zero
        inform%hard_case = .TRUE.
        DO i = 1, n
          IF ( H( i ) == lambda_min ) THEN
            IF ( ABS( C( i ) ) > epsmch * c_norm ) THEN
              inform%hard_case = .FALSE.
              c2 = c2 + C( i ) ** 2
            ELSE
              i_hard = i
            END IF
          END IF
        END DO

!  the hard case may occur

        IF ( inform%hard_case ) THEN
          DO i = 1, n
            IF ( H( i ) /= lambda_min ) THEN
              X( i )  = - C( i ) / ( H( i ) + lambda )
            ELSE
              X( i ) = zero
            END IF
          END DO
          inform%x_norm = TWO_NORM( X )

!  compute the target value ( lambda / sigma )^(1/(p-2))

          target = ( lambda / sigma ) ** oopm2

!  the hard case does occur

          IF ( inform%x_norm <= target ) THEN
            IF ( inform%x_norm < target ) THEN

!  compute the step alpha so that X + alpha E_i_hard lies on the trust-region
!  boundary and gives the smaller value of q

              utx = X( i_hard ) / target
              distx = ( target - inform%x_norm ) *                             &
                ( ( target + inform%x_norm ) / target )
              alpha = sign( distx / ( abs( utx ) +                             &
                            sqrt( utx ** 2 + distx / target ) ), utx )

!  record the optimal values

              X( i_hard ) = X( i_hard ) + alpha
            END IF
            inform%x_norm = TWO_NORM( X )
            inform%obj =                                                       &
                f + half * ( DOT_PRODUCT( C, X ) - lambda * target ** 2 )
            inform%obj_regularized = inform%obj + ( lambda / p ) * target ** 2

            IF ( printi ) THEN
              WRITE( out, "( A, A2, I4, 1X, 3ES22.15 )" ) prefix, region,      &
                it, inform%x_norm - target, lambda, ABS( delta_lambda )
              WRITE( out, "( A,                                                &
          &    ' Normal stopping criteria satisfied' )" ) prefix
            END IF
            inform%status = RAL_NLLS_ok
            GO TO 900

!  the hard case didn't occur after all

          ELSE
            inform%hard_case = .FALSE.

!  compute the first derivative of ||x|(lambda)||^2  ...

            w_norm2 = zero
            DO i = 1, n
              IF ( H( i ) /= lambda_min )                                      &
                w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
            END DO
            x_norm2( 1 ) = - two * w_norm2

!  ... and ( lambda / sigma )^(2/(p-2))

            IF ( p == three ) THEN
              theta_beta( 1 ) = two * lambda * oos2
            ELSE
              theta_beta( 1 ) =                                                &
                topm2 * ( lambda ** ( topm2 - one ) ) / ( sigma ** topm2 )
            END IF

!  compute the Newton correction

            lambda = lambda + ( inform%x_norm ** 2 - target ** 2 ) /           &
                              ( x_norm2( 1 ) - theta_beta( 1 ) )
            lambda_l = MAX( lambda_l, lambda )
          END IF

!  there is a singularity at lambda. Compute the point for which the
!  sum of squares of the singular terms is equal to target^2

        ELSE
          lambda = MAX( lambda * ( one + epsmch ),                             &
            DRQS_lambda_root( lambda_min, SQRT( c2 ) * sigma ** oopm2, oopm2 ) )
          lambda_l = MAX( lambda_l, lambda )
        END IF

!  lambda lies above the largest singularity.

      ELSE

!  compute the value of ||x(lambda)||

!       IF ( p == three ) THEN   !! For the time being, only p == 3
        IF ( .FALSE. ) THEN
          w_norm2 = zero
          DO i = 1, n
            w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 2
          END DO
          w_norm2 = SQRT( w_norm2 )

!  an upper bound on the required lambda occurs when this value is equal to
!  ( lambda / sigma )^(1/(p-2))

          ii = MINLOC( H )
          j = ii( 1)
          WRITE( out , * ) ' upper lambda = ', sigma * w_norm2, lambda_u
          lambda_u = MIN( sigma * w_norm2, lambda_u )

!  the function ||x(lambda)|| is no smaller than h(lambda)
!  c_j^2 / (lambda + lambda_j)^2 + sum_i/=j c_i^2 / (lambda_u + lambda_i)^2,
!  so the required lambda is no smaller than the largest root of
!  h(lambda) = ( lambda / sigma )^(2/(p-2))

          w_norm2 = zero
          DO i = 1, n
            IF ( i /= j )                                                      &
              w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda_u ) ** 2
          END DO
          w_norm2 = w_norm2 * sigma ** 2

          a_0 = - ( sigma * C( j ) ) ** 2 - w_norm2 * H( j ) ** 2
          a_1 =  - two * w_norm2 * H( j )
          a_2 = H( j ) ** 2 - w_norm2
          a_3 = two * H( j )
          a_4 = one
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ),                     &
                       ABS( a_3 ), ABS( a_4 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max
            a_2 = a_2 / a_max ; a_3 = a_3 / a_max ; a_4 = a_4 / a_max
          END IF
          CALL ROOTS_quartic( a_0, a_1, a_2, a_3, a_4, roots_tol, nroots,      &
                            roots( 1 ), roots( 2 ), roots( 3 ), roots( 4 ),    &
                            roots_debug )
          WRITE( out, * ) ' starting lambda = ', roots( : nroots )
!         lambda = MAX( roots( nroots ), lambda + epsmch )
          lambda = roots( nroots )
        END IF
      END IF

!  the iterates will all be in the L region. Prepare for the main loop

      it = 0
      max_order = MAX( 1, MIN( max_degree, control%taylor_max_degree ) )

!  start the main loop

!write(6,*) ' H ', H
!write(6,*) ' C ', C
!write(6,*) ' sigma ', sigma

      DO
        it = it + 1
        IF ( control%maxit > 0 .AND. it > control%maxit ) THEN
          inform%status = RAL_NLLS_error_max_iterations
          IF ( printi ) WRITE( out, "( ' iteration limit exceeded' )" )
          EXIT
        END IF

!  if H(lambda) is positive definite, solve  H(lambda) x = - c

        DO i = 1, n
          X( i )  = - C( i ) / ( H( i ) + lambda )
        END DO

!  compute the M-norm of x, ||x||_M

        inform%x_norm = TWO_NORM( X )
        x_norm2( 0 ) = inform%x_norm ** 2

!  compute the target value ( lambda / sigma )^(1/(p-2))

        target = ( lambda / sigma ) ** oopm2

!  the current estimate gives a good approximation to the required root

        IF ( ABS( inform%x_norm - target ) <=                                  &
             MAX( control%stop_normal * MAX( one, inform%x_norm, target ),     &
                  control%stop_absolute_normal  ) ) THEN
          IF ( inform%x_norm > target ) THEN
            region = 'L'
            lambda_l = MAX( lambda_l, lambda )
          ELSE
            region = 'G'
            lambda_u = MIN( lambda_u, lambda )
          END IF
          IF ( printt .AND. it > 1 ) WRITE( out, 2030 ) prefix
          IF ( printi ) THEN
            WRITE( out, "( A, A2, I4, 1X, 3ES22.15 )" ) prefix, region,        &
              it, inform%x_norm - target, lambda, ABS( delta_lambda )
            WRITE( out, "( A,                                                  &
        &    ' Normal stopping criteria satisfied' )" ) prefix
          END IF
          inform%status = RAL_NLLS_ok
          EXIT
        END IF

        lambda_l = MAX( lambda_l, lambda )

!  a lambda in L has been found. It is now simply a matter of applying
!  a variety of Taylor-series-based methods starting from this lambda

        IF ( printi ) WRITE( out, "( A, A2, I4, 1X, 3ES22.15 )" ) prefix,      &
          region, it, inform%x_norm - target, lambda, ABS( delta_lambda )

!  precaution against rounding producing lambda outside L

        IF ( lambda > lambda_u ) THEN
          inform%status = RAL_NLLS_error_ill_conditioned
          IF ( printi ) THEN
            WRITE( out, 2030 ) prefix
            WRITE( out, "( A, 2X, I4, 3ES22.15, /, A,                          &
           &               ' normal exit with lambda outside L' )" )           &
              prefix, it, inform%x_norm - target, lambda, ABS( delta_lambda ), &
              prefix
          END IF
          EXIT
        END IF

!  compute first derivatives of x^T M x

!  form ||w||^2 = x^T H^-1(lambda) x

        w_norm2 = zero
        DO i = 1, n
          w_norm2 = w_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 3
        END DO

!  compute the first derivative of x_norm2 = x^T M x

        x_norm2( 1 ) = - two * w_norm2

!  count the number of corrections computed

        n_lambda = 0

!  compute Taylor approximants of degree one;
!  special (but frequent) case when p = 3

        IF ( p == three ) THEN

!  compute pi_beta = ||x||^beta and its first derivative when beta = 2

          beta = two
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )

!  compute the Newton correction (for beta = 2)

          a_0 = pi_beta( 0 ) - target ** 2
          a_1 = pi_beta( 1 ) - two * lambda * oos2
          a_2 = - oos2
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max ; a_2 = a_2 / a_max
          END IF
          CALL ROOTS_quadratic( a_0, a_1, a_2, roots_tol, nroots,              &
                                roots( 1 ), roots( 2 ), roots_debug )
          lambda_plus = lambda + roots( 2 )
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF

!  compute pi_beta = ||x||^beta and its first derivative when beta = 1

          beta = one
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )

!  compute the Newton correction (for beta = 1)

          delta_lambda = - ( pi_beta( 0 ) - target ) / ( pi_beta( 1 ) - oos )
          lambda_plus = lambda + delta_lambda
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF

!  compute pi_beta = ||x||^beta and its first derivative when beta = - 1

          beta = - one
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )

!  compute the Newton correction (for beta = -1)

          a_0 = pi_beta( 0 ) * lambda - sigma
          a_1 = pi_beta( 0 ) + lambda * pi_beta( 1 )
          a_2 = pi_beta( 1 )
          a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ) )
          IF ( a_max > zero ) THEN
            a_0 = a_0 / a_max ; a_1 = a_1 / a_max ; a_2 = a_2 / a_max
          END IF
          CALL ROOTS_quadratic( a_0, a_1, a_2, roots_tol, nroots,              &
                                roots( 1 ), roots( 2 ), roots_debug )
          lambda_plus = lambda + roots( 2 )
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF

!  more general p

        ELSE

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their first derivatives when beta = p-2

          beta = pm2
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )
          CALL DRQS_theta_derivs( 1, beta / pm2, lambda, sigma,                &
                                  theta_beta( : 1 )  )

!  compute the "linear Taylor approximation" correction (for beta = p-2)

          delta_lambda = - ( pi_beta( 0 ) - theta_beta( 0 ) ) /                &
                           ( pi_beta( 1 ) - theta_beta( 1 ) )
          lambda_plus = lambda + delta_lambda
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their first derivatives when beta = (p-2)/2

          beta = pm2 / two
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )
          CALL DRQS_theta_derivs( 1, beta / pm2, lambda, sigma,                &
                                 theta_beta( : 1 )  )

!  compute the "linear Taylor approximation" correction (for beta = (p-2)/2)

          delta_lambda = - ( pi_beta( 0 ) - theta_beta( 0 ) ) /                &
                           ( pi_beta( 1 ) - theta_beta( 1 ) )
          lambda_plus = lambda + delta_lambda
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their first derivatives when beta = max(2-p,-1)

          beta = max( - pm2, - one )
          CALL DRQS_pi_derivs( 1, beta, x_norm2( : 1 ), pi_beta( : 1 ) )
          CALL DRQS_theta_derivs( 1, beta / pm2, lambda, sigma,                &
                                 theta_beta( : 1 ) )

!  compute the "linear Taylor approximation" correction (for beta = max(2-p,-1))

          delta_lambda = - ( pi_beta( 0 ) - theta_beta( 0 ) ) /                &
                           ( pi_beta( 1 ) - theta_beta( 1 ) )
          lambda_plus = lambda + delta_lambda
          IF (  lambda_plus < lambda ) THEN
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda_plus
          END IF
        END IF

        IF ( max_order >= 3 ) THEN

!  compute the second derivative of x^T x

          z_norm2 = zero
          DO i = 1, n
            z_norm2 = z_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 4
          END DO
          x_norm2( 2 ) = six * z_norm2

!  compute the third derivatives of x^T x

          v_norm2 = zero
          DO i = 1, n
            v_norm2 = v_norm2 + C( i ) ** 2 / ( H( i ) + lambda ) ** 5
          END DO
          x_norm2( 3 ) = - twentyfour * v_norm2

!  compute pi_beta = ||x||^beta and its derivatives for various beta
!  and the resulting Taylor series approximants

!  special (but frequent) case when p = 3

          IF ( p == three ) THEN

!  compute pi_beta = ||x||^beta and its derivatives when beta = 2

            beta = two
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )

!  compute the "cubic Taylor approximaton" step (beta = 2)

            a_0 = pi_beta( 0 ) - target ** 2
            a_1 = pi_beta( 1 ) - two * lambda * oos2
            a_2 = half * pi_beta( 2 ) - oos2
            a_3 = sixth * pi_beta( 3 )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                              roots( 1 ), roots( 2 ), roots( 3 ),              &
                              roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda + roots( 1 )

!  compute pi_beta = ||x||^beta and its derivatives when beta = 1

            beta = one
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )

!  compute the "cubic Taylor approximaton" step (beta = 1)

            a_0 = pi_beta( 0 ) - target
            a_1 = pi_beta( 1 ) - oos
            a_2 = half * pi_beta( 2 )
            a_3 = sixth * pi_beta( 3 )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                              roots( 1 ), roots( 2 ), roots( 3 ),              &
                              roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda + roots( 1 )

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^beta and
!  their derivatives when beta = - 0.4

            beta = - point4
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )
            CALL DRQS_theta_derivs( 3, beta, lambda, sigma,                    &
                                    theta_beta( : 3 )  )

!  compute the "cubic Taylor approximaton" step (beta = - 0.4)

            a_0 = pi_beta( 0 ) - theta_beta( 0 )
            a_1 = pi_beta( 1 ) - theta_beta( 1 )
            a_2 = half * ( pi_beta( 2 ) - theta_beta( 2 ) )
            a_3 = sixth * ( pi_beta( 3 ) - theta_beta( 3 ) )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                               roots( 1 ), roots( 2 ), roots( 3 ),             &
                               roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda + roots( 1 )

!  more general p

          ELSE

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their derivatives when beta = p-2

            beta = pm2
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )
            CALL DRQS_theta_derivs( 3, beta / pm2, lambda, sigma,              &
                                    theta_beta( : 3 )  )

!  compute the "cubic Taylor approximation" correction (for beta = p-2)

            a_0 = pi_beta( 0 ) - theta_beta( 0 )
            a_1 = pi_beta( 1 ) - theta_beta( 1 )
            a_2 = half * ( pi_beta( 2 ) - theta_beta( 2 ) )
            a_3 = sixth * ( pi_beta( 3 ) - theta_beta( 3 ) )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                              roots( 1 ), roots( 2 ), roots( 3 ),              &
                              roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda +                                  &
              DRQS_required_root( .TRUE., nroots, roots( : 3 ) )

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their derivatives when beta = (p-2)/2

            beta = pm2 / two
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )
            CALL DRQS_theta_derivs( 3, beta / pm2, lambda, sigma,              &
                                    theta_beta( : 3 )  )

!  compute the "cubic Taylor approximation" correction (for beta = (p-2)/2)

            a_0 = pi_beta( 0 ) - theta_beta( 0 )
            a_1 = pi_beta( 1 ) - theta_beta( 1 )
            a_2 = half * ( pi_beta( 2 ) - theta_beta( 2 ) )
            a_3 = sixth * ( pi_beta( 3 ) - theta_beta( 3 ) )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                              roots( 1 ), roots( 2 ), roots( 3 ),              &
                              roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda +                                  &
              DRQS_required_root( .TRUE., nroots, roots( : 3 ) )

!  compute pi_beta = ||x||^beta and theta_beta = (lambda/sigma)^(beta/(p-2)) and
!  their derivatives when beta = max(2-p,-0.4)

            beta = max( - pm2, - point4 )
            CALL DRQS_pi_derivs( 3, beta, x_norm2( : 3 ), pi_beta( : 3 ) )
            CALL DRQS_theta_derivs( 3, beta / pm2, lambda, sigma,              &
                                    theta_beta( : 3 )  )

!  compute the "cubic Taylor approximation" correction (for beta=max(2-p,-0.4))

            a_0 = pi_beta( 0 ) - theta_beta( 0 )
            a_1 = pi_beta( 1 ) - theta_beta( 1 )
            a_2 = half * ( pi_beta( 2 ) - theta_beta( 2 ) )
            a_3 = sixth * ( pi_beta( 3 ) - theta_beta( 3 ) )
            a_max = MAX( ABS( a_0 ), ABS( a_1 ), ABS( a_2 ), ABS( a_3 ) )
            IF ( a_max > zero ) THEN
              a_0 = a_0 / a_max ; a_1 = a_1 / a_max
              a_2 = a_2 / a_max ; a_3 = a_3 / a_max
            END IF
            CALL ROOTS_cubic( a_0, a_1, a_2, a_3, roots_tol, nroots,           &
                              roots( 1 ), roots( 2 ), roots( 3 ),              &
                              roots_debug )
            n_lambda = n_lambda + 1
            lambda_new( n_lambda ) = lambda +                                  &
              DRQS_required_root( .TRUE., nroots, roots( : 3 ) )
          END IF
        END IF

!  record all of the estimates of the optimal lambda

        IF ( printd ) THEN
          MLOC(:) = MAXLOC( lambda_new( : n_lambda ) )
          WRITE( out, "( A, ' lambda_t (', I1, ')', 3ES20.13 )" )              &
            prefix, MLOC(1), lambda_new( : MIN( 3, n_lambda ) )
          IF ( n_lambda > 3 ) WRITE( out, "( A, 13X, 3ES20.13 )" )             &
            prefix, lambda_new( 4 : MIN( 6, n_lambda ) )
          IF ( n_lambda > 6 ) WRITE( out, "( A, 13X, 3ES20.13 )" )             &
            prefix, lambda_new( 7 : MIN( 9, n_lambda ) )
        END IF

!  compute the best Taylor improvement

        lambda_plus = MAXVAL( lambda_new( : n_lambda ) )
        delta_lambda = lambda_plus - lambda
        lambda = lambda_plus

!  improve the lower bound if possible

        lambda_l = MAX( lambda_l, lambda_plus )

!  check that the best Taylor improvement is significant

!write(6,*) ABS( delta_lambda ), epsmch * MAX( one, ABS( lambda ) )
!       IF ( ABS( delta_lambda ) < epsmch * MAX( one, ABS( lambda ) ) ) THEN
        IF ( ABS( delta_lambda ) < epsmch * ABS( lambda ) ) THEN
          inform%status = RAL_NLLS_ok
          IF ( printi ) WRITE( out, "( A, ' normal exit with no ',             &
         &                     'significant Taylor improvement' )" ) prefix
          EXIT
        END IF

!  End of main iteration loop

      END DO

!  Record the optimal obective value

      inform%obj = f + half * ( DOT_PRODUCT( C, X ) - lambda * target ** 2 )
      inform%obj_regularized = inform%obj + ( lambda / p ) * target ** 2
      IF ( printi ) WRITE( out,                                                &
        "( A, ' estimated, true objective values =', 2ES21.13 )" ) prefix,     &
          inform%obj_regularized, f + DOT_PRODUCT( C, X ) +                    &
            half * DOT_PRODUCT( X, H( : n ) * X ) +                            &
          ( sigma / p ) * inform%x_norm ** p

!  ----
!  Exit
!  ----

 900  CONTINUE
      inform%multiplier = lambda
      inform%pole = MAX( zero, - lambda_min )
      RETURN


! Non-executable statements

 2030 FORMAT( A, '    it     ||x||-target              lambda ',               &
                 '              d_lambda' )
!2050 FORMAT( A, ' time( SLS_solve ) = ', F0.2 )

!  End of subroutine DRQS_solve_main

      END SUBROUTINE DRQS_solve_main

!-*-*-*-*-*-  D R Q S _ P I _ D E R I V S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE DRQS_pi_derivs( max_order, beta, x_norm2, pi_beta )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute pi_beta = ||x||^beta and its derivatives
!
!  Arguments:
!  =========
!
!  Input -
!   max_order - maximum order of derivative
!   beta - power
!   x_norm2 - (0) value of ||x||^2,
!             (i) ith derivative of ||x||^2, i = 1, max_order
!  Output -
!   pi_beta - (0) value of ||x||^beta,
!             (i) ith derivative of ||x||^beta, i = 1, max_order
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: max_order
      REAL ( KIND = wp ), INTENT( IN ) :: beta, x_norm2( 0 : max_order )
      REAL ( KIND = wp ), INTENT( OUT ) :: pi_beta( 0 : max_order )

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      REAL ( KIND = wp ) :: hbeta

      hbeta = half * beta
      pi_beta( 0 ) = x_norm2( 0 ) ** hbeta
      IF ( hbeta == one ) THEN
        pi_beta( 1 ) = x_norm2( 1 )
        IF ( max_order == 1 ) RETURN
        pi_beta( 2 ) = x_norm2( 2 )
        IF ( max_order == 2 ) RETURN
        pi_beta( 3 ) = x_norm2( 3 )
      ELSE IF ( hbeta == two ) THEN
        pi_beta( 1 ) = two * x_norm2( 0 ) * x_norm2( 1 )
        IF ( max_order == 1 ) RETURN
        pi_beta( 2 ) = two * ( x_norm2( 1 ) ** 2 + x_norm2( 0 ) * x_norm2( 2 ) )
        IF ( max_order == 2 ) RETURN
        pi_beta( 3 ) = two *                                                   &
          ( x_norm2( 0 ) * x_norm2( 3 ) + three * x_norm2( 1 ) * x_norm2( 2 ) )
      ELSE
        pi_beta( 1 )                                                           &
          = hbeta * ( x_norm2( 0 ) ** ( hbeta - one ) ) * x_norm2( 1 )
        IF ( max_order == 1 ) RETURN
        pi_beta( 2 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - two ) ) *           &
          ( ( hbeta - one ) * x_norm2( 1 ) ** 2 + x_norm2( 0 ) * x_norm2( 2 ) )
        IF ( max_order == 2 ) RETURN
        pi_beta( 3 ) = hbeta * ( x_norm2( 0 ) ** ( hbeta - three ) ) *         &
          ( x_norm2( 3 ) * x_norm2( 0 ) ** 2 + ( hbeta - one ) *               &
            ( three * x_norm2( 0 ) * x_norm2( 1 ) * x_norm2( 2 ) +             &
              ( hbeta - two ) * x_norm2( 1 ) ** 3 ) )
      END IF
      RETURN

!  End of subroutine DRQS_pi_derivs

      END SUBROUTINE DRQS_pi_derivs

!-*-*-*-*-  D R Q S _ T H E T A _ D E R I V S   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE DRQS_theta_derivs( max_order, beta, lambda, sigma, theta_beta )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute theta_beta = (lambda/sigma)^beta and its derivatives
!
!  Arguments:
!  =========
!
!  Input -
!   max_order - maximum order of derivative
!   beta - power
!   lambda, sigma - lambda and sigma
!  Output -
!   theta_beta - (0) value of (lambda/sigma)^beta,
!             (i) ith derivative of (lambda/sigma)^beta, i = 1, max_order
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: max_order
      REAL ( KIND = wp ), INTENT( IN ) :: beta, lambda, sigma
      REAL ( KIND = wp ), INTENT( OUT ) :: theta_beta( 0 : max_order )

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      REAL ( KIND = wp ) :: los, oos

      los = lambda / sigma
      oos = one / sigma

      theta_beta( 0 ) = los ** beta
      IF ( beta == one ) THEN
        theta_beta( 1 ) = oos
        IF ( max_order == 1 ) RETURN
        theta_beta( 2 ) = zero
        IF ( max_order == 2 ) RETURN
        theta_beta( 3 ) = zero
      ELSE IF ( beta == two ) THEN
        theta_beta( 1 ) = two * los * oos
        IF ( max_order == 1 ) RETURN
        theta_beta( 2 ) = oos ** 2
        IF ( max_order == 2 ) RETURN
        theta_beta( 3 ) = zero
      ELSE
        theta_beta( 1 ) = beta * ( los ** ( beta - one ) ) * oos
        IF ( max_order == 1 ) RETURN
        theta_beta( 2 ) = beta * ( los ** ( beta - two ) ) *                   &
                          ( beta - one ) * oos ** 2
        IF ( max_order == 2 ) RETURN
        theta_beta( 3 ) = beta * ( los ** ( beta - three ) ) *                 &
                          ( beta - one ) * ( beta - two ) * oos ** 3
      END IF

      RETURN

!  End of subroutine DRQS_theta_derivs

      END SUBROUTINE DRQS_theta_derivs

!-*-*-*-*-*-  D R Q S _ R E Q U I R E D _ R O O T  F U C T I O N   -*-*-*-*-*-

      FUNCTION DRQS_required_root( positive, nroots, roots )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Determine the required root of the three roots of the secular equation.
!  This is either the most positive root (positive=.TRUE.) or the least
!  negative one (positive=.FALSE.)
!
!  Arguments:
!  =========
!
!  Input -
!   positive - .TRUE. if the largest positive root is required,
!               .FALSE. if the least negative one
!   nroots - number of roots
!   roots - roots in increasing order
!  Output -
!   DRQS_required root - the required root
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      REAL ( KIND = wp ) :: DRQS_required_root
      INTEGER, INTENT( IN ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: roots
      LOGICAL, INTENT( IN ) :: positive

      IF ( positive ) THEN
        IF ( SIZE( roots ) == 3 ) THEN
          IF ( nroots == 3 ) THEN
            DRQS_required_root = roots( 3 )
          ELSE IF ( nroots == 2 ) THEN
            DRQS_required_root = roots( 2 )
          ELSE
            DRQS_required_root = roots( 1 )
          END IF
        ELSE
          IF ( nroots == 2 ) THEN
            DRQS_required_root = roots( 2 )
          ELSE
            DRQS_required_root = roots( 1 )
          END IF
        END IF
      ELSE
        IF ( SIZE( roots ) == 3 ) THEN
          IF ( nroots == 3 ) THEN
            IF ( roots( 3 ) > zero ) THEN
              IF ( roots( 2 ) > zero ) THEN
                DRQS_required_root = roots( 1 )
              ELSE
                DRQS_required_root = roots( 2 )
              END IF
            ELSE
              DRQS_required_root = roots( 3 )
            END IF
          ELSE IF ( nroots == 2 ) THEN
            IF ( roots( 2 ) > zero ) THEN
              DRQS_required_root = roots( 1 )
            ELSE
              DRQS_required_root = roots( 2 )
            END IF
          ELSE
            DRQS_required_root = roots( 1 )
          END IF
        ELSE
          IF ( nroots == 2 ) THEN
            IF ( roots( 2 ) > zero ) THEN
              DRQS_required_root = roots( 1 )
            ELSE
              DRQS_required_root = roots( 2 )
            END IF
          ELSE
            DRQS_required_root = roots( 1 )
          END IF
        END IF
      END IF
      RETURN

!  End of function DRQS_required_root

      END FUNCTION DRQS_required_root

!-*-*-*-*-*-*-  D R Q S _ L A M B D A  _ R O O T  F U C T I O N   -*-*-*-*-*-*-

      FUNCTION DRQS_lambda_root( a, b, power )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the positive root of lambda + a = b/lambda^power
!
!  Arguments:
!  =========
!
!  Input -
!   a, b, power - data for the above problem (with b, power > 0)
!  Output -
!   DRQS_lambda root - the required root
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      REAL ( KIND = wp ) :: DRQS_lambda_root
      REAL ( KIND = wp ), INTENT( IN ) :: a, b, power

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      INTEGER :: nroots, it
!     INTEGER, PARAMETER :: newton_max = 10000
      INTEGER, PARAMETER :: newton_max = 20
      REAL ( KIND = wp ) :: lambda, phi, phip, d_lambda, other, power_plus_1

!write(6,*) ' a, b, p', a, b, power

!  special case: a = 0 = b

      IF ( a == zero .AND. b == zero ) THEN
        DRQS_lambda_root = zero ; RETURN
      END IF

      power_plus_1 = power + one

!  compute as initial lower bound on the root

      IF ( power == one ) THEN
        CALL ROOTS_quadratic( - b , a, one, roots_tol, nroots, other, lambda,  &
                              roots_debug )
      ELSE

!  when power > 1, 1/lambda <= 1/lambda^p for lambda in (0,1]

        IF ( power > one ) THEN
          CALL ROOTS_quadratic( - b , a, one, roots_tol, nroots, other,        &
                                lambda, roots_debug )
          lambda = MIN( one, lambda )
        ELSE
          lambda = epsmch
        END IF

!  check if lambda = 1 is acceptable

        IF ( one + a <= b ) lambda = MAX( lambda, one )

!  when a > 0, find where the tangent to b/lambda^power at
!  lambda = b^(1/power+1) intersects lambda + a

        IF ( a >= zero ) THEN
          lambda = MAX( lambda, b ** ( one / power_plus_1 ) - a / power_plus_1 )

!  when a < 0, both the lambda-intercept of lambda + a and the interection
!  of lambda with beta / lambda^(1/power+1) give lower bounds on the root

        ELSE
          lambda = MAX( lambda, - a, b ** ( one / power_plus_1 ) )
        END IF

!  perform Newton's method to refine the root

        DO it = 1, newton_max
          phi = lambda + a - b / ( lambda ** power )
!         write(6,*) ' lambda ', lambda, phi
          IF ( ABS( phi ) <= ten * epsmch *                                   &
                 MAX(  lambda + a, b / ( lambda ** power ) ) ) EXIT
          phip = one + b * power / ( lambda ** power_plus_1 )
          d_lambda = - phi / phip
          IF ( ABS( d_lambda ) <= epsmch * MAX( one, lambda ) ) EXIT
          lambda = lambda + d_lambda
        END DO
      END IF
      DRQS_lambda_root = lambda

      RETURN

!  End of function DRQS_lambda_root

      END FUNCTION DRQS_lambda_root

!-*-*-*-*-*-*-*-*-*-  T W O  _ N O R M   F U N C T I O N   -*-*-*-*-*-*-*-*-*-

       FUNCTION TWO_NORM( X )

!  Compute the l_2 norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: TWO_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         TWO_NORM = NRM2( n, X, 1 )
       ELSE
         TWO_NORM = zero
       END IF
       RETURN

!  End of function TWO_NORM

       END FUNCTION TWO_NORM

!-*-*-*-*-*-  End of R A L _ N L L S _ R Q S  double  M O D U L E  *-*-*-*-*-*-

   END MODULE MODULE_PREC(RAL_NLLS_DRQS)
