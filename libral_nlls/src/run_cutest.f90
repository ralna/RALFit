program run_cutest

  use iso_c_binding
  use nlls_module

  ! A local (and more customizable) version of ral_nlls_main
  ! in order to run the test problems....
  !
  ! Tyrone Rees, 2015 (heavily based on Nick Gould's version in cutest)

  implicit none

  type, extends( params_base_type ) :: user_type
     ! still empty
end type user_type

INTEGER :: status, m, n
INTEGER( c_int ) :: len_work_integer, len_work_real
REAL( c_double ), PARAMETER :: infty = 1.0D+19
REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: X, X_l, X_u
REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: Y, C_l, C_u, F
REAL( c_double ), DIMENSION( : ), ALLOCATABLE ::  Work_real
INTEGER( c_int ), DIMENSION( : ), ALLOCATABLE ::  Work_integer
type( user_type ), target :: params
TYPE( NLLS_inform_type ) :: inform
TYPE( NLLS_control_type ) :: control
LOGICAL, DIMENSION( : ), ALLOCATABLE  :: EQUATN, LINEAR
CHARACTER ( LEN = 10 ) :: pname
CHARACTER ( LEN = 30 ) :: progressfile
CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: VNAMES, CNAMES
REAL( c_double ), DIMENSION( 2 ) :: CPU
REAL( c_double ), DIMENSION( 4 ) :: CALLS
INTEGER :: io_buffer = 11
INTEGER, PARAMETER :: input = 55, indr = 46, out = 6
integer :: i
logical :: read_controls_in
!  Interface blocks

INTERFACE
  SUBROUTINE eval_F( status, n, m, X, F, params )
    USE ISO_C_BINDING
    import :: params_base_type
    INTEGER ( c_int ), INTENT( OUT ) :: status
    INTEGER ( c_int ), INTENT( IN ) :: n, m
    REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
    REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: F
    class( params_base_type ), intent(in) :: params
  END SUBROUTINE eval_F
END INTERFACE

INTERFACE
  SUBROUTINE eval_J( status, n, m, X, J, params )
    USE ISO_C_BINDING
    import :: params_base_type
    INTEGER ( c_int ), INTENT( OUT ) :: status
    INTEGER ( c_int ), INTENT( IN ) :: n, m
    REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
    REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: J
    class( params_base_type ), intent(in) :: params
  END SUBROUTINE eval_J
END INTERFACE

INTERFACE
  SUBROUTINE eval_HF( status, n, m, X, F, H, params )
    USE ISO_C_BINDING
    import :: params_base_type
    INTEGER ( c_int ), INTENT( OUT ) :: status
    INTEGER ( c_int ), INTENT( IN ) :: n, m
    REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
    REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: F
    REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: H
    class( params_base_type ), intent(in) :: params
  END SUBROUTINE eval_HF
END INTERFACE

!  open the relevant file

OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
REWIND( input )

!  compute problem dimensions

CALL CUTEST_cdimen( status, input, n, m )
IF ( status /= 0 ) GO TO 910

!  allocate space 

ALLOCATE( X( n ), X_l( n ), X_u( n ), Y( m ), C_l( m ), C_u( m ),        &
    EQUATN( m ), LINEAR( m ), STAT = status )
IF ( status /= 0 ) GO TO 990

!  initialize problem data structure

!  set up the data structures necessary to hold the problem functions.

CALL CUTEST_csetup( status, input, out, io_buffer, n, m,                 &
    X, X_l, X_u, Y, C_l, C_u, EQUATN, LINEAR, 0, 0, 0 )
IF ( status /= 0 ) GO TO 910
CLOSE( input )

!  allocate more space 

DEALLOCATE( X_l, X_u, Y, C_l, C_u, EQUATN, LINEAR )
len_work_integer = 0
len_work_real = m + n * ( m + n )
ALLOCATE( Work_integer( len_work_integer ), Work_real( len_work_real ),  &
    STAT = status )
IF ( status /= 0 ) GO TO 990

!  open the Spec file for the method

OPEN( indr, FILE = '../cutest/sif/RAL_NLLS.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
REWIND( indr )

!  read input Spec data

!  error = unit for error messages
!  out = unit for information

!  set up algorithmic input data

READ ( indr, 1000 ) control%error, control%out, control%print_level
CLOSE ( indr )

!  call the minimizer
read_controls_in = .true.
if (read_controls_in) then
   open(unit = 17, file = "cutest_control.in")
   read(17,*) control%print_level
   read(17,*) control%nlls_method
   read(17,*) control%subproblem_eig_fact
   read(17,*) control%maxit
   read(17,*) control%model
   read(17,*) control%stop_g_relative, control%stop_g_absolute
   read(17,*) control%output_progress_vectors
   close(unit=17)
else
   control%print_level = 3
   control%nlls_method = 3
   control%subproblem_eig_fact = .true.
   control%maxit = 1000
   control%model = 7
   control%stop_g_relative = 1e-12
   control%stop_g_absolute = 1e-12
   control%output_progress_vectors = .true.
end if

inform%iter = 23
open(unit=42,file="data/results.out",position="append")
open(unit=52,file="data/results.out_data",position="append")
CALL RAL_NLLS( n, m, X,                                    &
     !              Work_integer, len_work_integer, Work_real, len_work_real &
    eval_F, eval_J, eval_HF, params,            &
    inform, control)!, inform )
IF ( status /= 0 ) GO TO 910

!  output report

CALL CUTEST_ureport( status, CALLS, CPU )
IF ( status /= 0 ) GO TO 910

ALLOCATE( F( m ), VNAMES( n ), CNAMES( m ), STAT = status )
CALL CUTEST_cnames( status, n, m, pname, VNAMES, CNAMES )

write(42,'(a,a,i0,a,i0,a,i0,a,i0,a,ES12.4)') pname,': n = ', n, ', m = ',m, &
     ',   status = ', inform%status, &                                   
     ',   iter = ', inform%iter , &
     ',   obj = ', inform%obj
write(52,'(i0,T10,i0,T20,i0,T30,i0,T40,ES23.15E3,T65,ES23.15E3)') n,&
                                                  m,&
                                                  inform%status, &
                                                  inform%iter, &
                                                  inform%obj, & 
                                                  inform%norm_g
write(*,'(a,a,i0,a,i0,a,i0,a,i0,a,ES12.4)') pname,': n = ', n, ', m = ',m, &
     ',   status = ', inform%status, &                                   
     ',   iter = ', inform%iter, &
     ',   obj = ', inform%obj
close(unit=42)


if (control%output_progress_vectors) then
   progressfile = "data/" // trim(pname) // "_progress.out"
   open(unit = 62, file = progressfile)
   do i = 1,inform%iter + 1
      write(62,'(ES20.14, T25, ES20.14)') inform%gradvec(i), &
                                          inform%resvec(i)
   end do
   close(unit = 62)
end if



CALL eval_F( status, n, m, X, F, params)

!WRITE( out, 2110 ) ( i, VNAMES( i ), X( i ), i = 1, n )
!WRITE( out, 2120 ) ( i, CNAMES( i ), F( i ), i = 1, m )
!!WRITE( out, 2000 ) pname, n, CALLS( 1 ), inform%obj, CPU( 1 ), CPU( 2 )
!WRITE( out, 2000 ) pname, n, CALLS( 1 ), 100.0, CPU( 1 ), CPU( 2 )

!  clean-up data structures

DEALLOCATE( X, F, VNAMES, CNAMES, Work_integer, Work_real,               &
    STAT = status )
IF ( status /= 0 ) GO TO 910
CALL CUTEST_cterminate( status )
STOP

!  error returns

910 CONTINUE
WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
STOP

990 CONTINUE
WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
STOP

!  Non-executable statements

!!$2000 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //,                    &
!!$    ' Package used            :  RAL_NLLS ',  /,                         &
!!$    ' Problem                 :  ', A10,    /,                           &
!!$    ' # variables             =      ', I10 /,                           &
!!$    ' # residuals             =        ', F8.2 /,                        &
!!$    ' Final f                 = ', E15.7 /,                              &
!!$    ' Set up time             =      ', 0P, F10.2, ' seconds' /,         &
!!$    ' Solve time              =      ', 0P, F10.2, ' seconds' //,        &
!!$    66('*') / )
1000 FORMAT( I6, /, I6, /, I6 )
!!$2110 FORMAT( /, ' The variables:', /, &
!!$    '     i name          value',  /, ( I6, 1X, A10, 1P, D12.4 ) )
!!$2120 FORMAT( /, ' The constraints:', /, '     i name          value',          &
!!$    /, ( I6, 1X, A10, 1P, D12.4 ) )

!  End of RAL_NLLS_main



end program


SUBROUTINE eval_F( status, n, m, X, F, params )
USE ISO_C_BINDING
use :: nlls_module, only : params_base_type

INTEGER ( c_int ), INTENT( OUT ) :: status
INTEGER ( c_int ), INTENT( IN ) :: n, m
REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
REAL ( c_double ), DIMENSION( m ), INTENT( OUT ) :: F
class( params_base_type ), intent(in) :: params
REAL ( c_double ) :: obj

!  evaluate the residuals F

CALL CUTEST_cfn( status, n, m, X, obj, F )
RETURN
END SUBROUTINE eval_F

SUBROUTINE eval_J( status, n, m, X, J, params)
USE ISO_C_BINDING
use :: nlls_module, only : params_base_type

INTEGER ( c_int ), INTENT( OUT ) :: status
INTEGER ( c_int ), INTENT( IN ) :: n, m
REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
REAL ( c_double ), DIMENSION( m * n ), INTENT( OUT ) :: J
class( params_base_type ), intent(in) :: params
REAL ( c_double ), DIMENSION( n ) :: G
REAL ( c_double ), DIMENSION( m ) :: Y
REAL ( c_double ), DIMENSION( m , n ) :: Jmatrix

!  evaluate the residual Jacobian J

CALL CUTEST_cgr( status, n, m, X, Y, .FALSE., G, .FALSE., m, n, Jmatrix ) 
! convert the Jacobian to a vector....
J = reshape(Jmatrix, (/n*m/) )
RETURN
END SUBROUTINE eval_J

SUBROUTINE eval_HF( status, n, m, X, F, H, params)
USE ISO_C_BINDING
use :: nlls_module, only : params_base_type

INTEGER ( c_int ), INTENT( OUT ) :: status
INTEGER ( c_int ), INTENT( IN ) :: n, m
REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
REAL ( c_double ), DIMENSION( m ), INTENT( IN ) :: F
REAL ( c_double ), DIMENSION( n*n ), INTENT( OUT ) :: H
class( params_base_type ), intent(in) :: params

real ( c_double ), dimension(n,n) :: Hmatrix
!  evaluate the product H = sum F_i Hessian F_i

CALL CUTEST_cdhc( status, n, m, X, F, n, Hmatrix )
H = reshape(Hmatrix, (/n*n/) )
RETURN
END SUBROUTINE eval_HF




