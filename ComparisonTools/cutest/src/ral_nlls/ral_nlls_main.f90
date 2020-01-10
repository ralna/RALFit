!   ( Last modified on 6 Oct 2015 at 15:50:00 )

      PROGRAM RAL_NLLS_main
      USE ISO_C_BINDING
      USE RAL_NLLS_DOUBLE

!  RAL_NLLS test driver for problems derived from SIF files

!  Nick Gould, October 2015

      IMPLICIT NONE

      type, extends( params_base_type ) :: user_type
         logical :: GotH
         integer :: lchp
         integer(c_int), allocatable :: chp_ptr(:)
         integer(c_int), allocatable :: chp_ind(:)
         real(c_double), allocatable :: chp_val(:)
      end type user_type

      INTEGER :: status, i, m, n
      REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: X, X_l, X_u
      REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: Y, C_l, C_u, F
      type( user_type ), target :: params
      TYPE( Nlls_inform ) :: inform
      TYPE( Nlls_options ) :: control
      LOGICAL, DIMENSION( : ), ALLOCATABLE  :: EQUATN, LINEAR
      CHARACTER ( LEN = 10 ) :: pname
      CHARACTER ( LEN = 30 ) :: summary_file = REPEAT( ' ', 30 )
      CHARACTER ( LEN = 30 ) :: iter_summary_file = REPEAT( ' ', 30 )
      CHARACTER ( LEN = 20 ) :: control_name = REPEAT( ' ', 20 ) 
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: VNAMES, CNAMES
      REAL( c_double ), DIMENSION( 2 ) :: CPU
      REAL( c_double ), DIMENSION( 7 ) :: CALLS
      INTEGER :: io_buffer = 11
      INTEGER :: summary_unit, iter_summary_unit, sol_unit, iores
      INTEGER, PARAMETER :: input = 55, indr = 46, out = 6
      LOGICAL :: filexx
      INTEGER :: fnevals, jacevals, hessevals, localiter, inner_iter
      REAL( c_double ) :: solve_time
      integer :: supply_eval_hp

!  open the relevant files

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

!  open the Spec file for the method

      OPEN( indr, FILE = 'RAL_NLLS.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
      REWIND( indr )

!  read input Spec data

!  error = unit for error messages
!  out = unit for information
!  print_level = level of output desired (<=0 = nothing)
!  nlls_method = method used (1=dogleg, 2=AINT, 3=More-Sorensen)
!  model = model used (1=first order, 2=Newton)
!  initial_radius = initial TR radius
!  stop_g_absolute = absolute stopping tolerance
!  stop_g_relative = relative stopping tolerance
!  summary_unit = write a one line summary to this unit (-ve = don't write)
!  summary_file = file name for summary (20 chars max)

!  set up algorithmic input data
      READ( indr, "(A)") control_name
      READ( indr, "(I6)") control%out
!      READ( indr, "(I6)") control%error
      READ( indr, "(I6)") control%print_level
      READ( indr, "(L)") control%print_options
      READ( indr, "(I6)") control%print_header
      READ( indr, "(I6)") control%maxit
      READ( indr, "(I6)") control%model
      READ( indr, "(I6)") control%type_of_method                                             
      READ( indr, "(I6)") control%nlls_method
      READ( indr, "(I6)") control%lls_solver
      READ( indr, "(E12.0)") control%stop_g_absolute
      READ( indr, "(E12.0)") control%stop_g_relative
      READ( indr, "(E12.0)") control%stop_f_absolute
      READ( indr, "(E12.0)") control%stop_f_relative
      READ( indr, "(E12.0)") control%stop_s
      READ( indr, "(I6)") control%relative_tr_radius
      READ( indr, "(E12.0)") control%initial_radius_scale
      READ( indr, "(E12.0)") control%initial_radius
      READ( indr, "(E12.0)") control%base_regularization
      READ( indr, "(I6)") control%regularization
      READ( indr, "(E12.0)") control%regularization_term
      READ( indr, "(E12.0)") control%regularization_power
      READ( indr, "(E12.0)") control%maximum_radius
      READ( indr, "(E12.0)") control%eta_successful
      READ( indr, "(E12.0)") control%eta_success_but_reduce
      READ( indr, "(E12.0)") control%eta_very_successful
      READ( indr, "(E12.0)") control%eta_too_successful
      READ( indr, "(E12.0)") control%radius_increase
      READ( indr, "(E12.0)") control%radius_reduce
      READ( indr, "(E12.0)") control%radius_reduce_max
      READ( indr, "(I6)") control%tr_update_strategy
      READ( indr, "(E12.0)") control%hybrid_switch
      READ( indr, "(L)") control%exact_second_derivatives
      READ( indr, "(L)") control%subproblem_eig_fact
      READ( indr, "(L)") control%use_ews_subproblem
      READ( indr, "(I6)") control%scale
      READ( indr, "(E12.0)") control%scale_max
      READ( indr, "(E12.0)") control%scale_min
      READ( indr, "(L)") control%scale_trim_min
      READ( indr, "(L)") control%scale_trim_max
      READ( indr, "(L)") control%scale_require_increase
      READ( indr, "(L)") control%setup_workspaces
      READ( indr, "(L)") control%remove_workspaces
      READ( indr, "(I6)") control%more_sorensen_maxits
      READ( indr, "(E12.0)") control%more_sorensen_shift
      READ( indr, "(E12.0)") control%more_sorensen_tiny
      READ( indr, "(E12.0)") control%more_sorensen_tol
      READ( indr, "(E12.0)") control%hybrid_tol
      READ( indr, "(I6)") control%hybrid_switch_its                                          
      READ( indr, "(E12.0)") control%reg_order
      READ( indr, "(I6)") control%inner_method
      READ( indr, "(L)") control%output_progress_vectors
      READ( indr, "(L)") control%update_lower_order
      READ( indr, "(I6)") summary_unit
      READ( indr, "(I6)") iter_summary_unit
      READ( indr, "(I6)") supply_eval_hp

!!$      READ( indr, "( I6, 8( /, I6 ), 11( /, E12.0 ), 3( /, I6 ) 2( /, E12.0),  &
!!$                     5 ( /, L20 ), /, I6, /, A, /, I6, /, A ) ")               &
!!$           control%error,                                                      &
!!$           control%out,                                                        &
!!$           control%print_level,                                                &
!!$           control%maxit,                                                      &
!!$           control%model,                                                      &
!!$           control%type_of_method,                                             &
!!$           control%nlls_method,                                                &
!!$           control%lls_solver,                                                 &
!!$           control%stop_g_absolute,                                            &
!!$           control%stop_g_relative,                                            &
!!$           control%stop_f_absolute,                                            &
!!$           control%stop_f_relative,                                            &
!!$           control%stop_s,                                                     &
!!$           control%relative_tr_radius,                                         &
!!$           control%initial_radius_scale,                                       &
!!$           control%initial_radius,                                             &
!!$           control%base_regularization,                                        &
!!$           control%regularized,                                                &
!!$           control%regularization_weight,                                      &
!!$           control%regularization_power,                                       &
!!$           control%maximum_radius,                                             &
!!$           control%eta_successful,                                             &
!!$           control%eta_success_but_reduce,                                     &
!!$           control%eta_very_successful,                                        &
!!$           control%eta_too_successful,                                         &
!!$           control%radius_increase,                                            &
!!$           control%radius_reduce,                                              &
!!$           control%radius_reduce_max,                                          &
!!$           control%tr_update_strategy,                                         &
!!$           control%hybrid_switch,                                              &
!!$           control%exact_second_derivatives,                                   &
!!$           control%subproblem_eig_fact,                                        &
!!$           control%scale,                                                      &
!!$           control%scale_min,                                                  &
!!$           control%scale_max,                                                  &
!!$           control%scale_trim_min,                                             &
!!$           control%scale_trim_max,                                             &
!!$           control%scale_require_increase,                                     &
!!$           control%calculate_svd_J,                                            &
!!$           control%setup_workspaces,                                           &
!!$           control%remove_workspaces,                                          &
!!$           control%more_sorensen_maxits,                                       &
!!$           control%more_sorensen_shift,                                        &
!!$           control%more_sorensen_tiny,                                         &
!!$           control%more_sorensen_tol,                                          &
!!$           control%hybrid_tol,                                                 &
!!$           control%hybrid_switch_its,                                          &
!!$           control%inner_method,                                               &
!!$           control%output_progress_vectors,                                    &
!!$           control%exact_second_derivatives,                                   &
!!$           summary_unit,                                                       &
!!$           summary_file,                                                       &
!!$           iter_summary_unit,                                                  &
!!$           iter_summary_file             

      CLOSE ( indr )
                       
      summary_file = trim(control_name)//".out" ! READ( indr, "(A)") summary_file      
      iter_summary_file = trim(control_name)//"_iter.out"
      write(*,*) 'summary_file = ', summary_file
      write(*,*) 'iter_summary_file = ', iter_summary_file
      write(*,*) 'control%print_level = ', control%print_level
      write(*,*) 'control%maxit = ', control%maxit
      write(*,*) 'control%scale = ', control%scale
      write(*,*) 'control%more_sorensen_maxits = ', control%more_sorensen_maxits
      write(*,*) 'control%output_progress_vectors = ', control%output_progress_vectors
      write(6,*) summary_unit, summary_file
      !write(6,*) iter_summary_unit, iter_summary_file
      write(6,*) iter_summary_unit, iter_summary_file


      IF ( summary_unit > 0 ) THEN
        INQUIRE( FILE = summary_file, EXIST = filexx )
        IF ( filexx ) THEN
           OPEN( summary_unit, FILE = summary_file, FORM = 'FORMATTED',        &
               STATUS = 'OLD', IOSTAT = iores , position="append")
        ELSE
           OPEN( summary_unit, FILE = summary_file, FORM = 'FORMATTED',        &
                STATUS = 'NEW', IOSTAT = iores )
        END IF
        IF ( iores /= 0 ) THEN 
          write( out, "( ' IOSTAT = ', I0, ' when opening file ', A,           &
        &  '. Stopping ' )" ) iores, summary_file
          STOP
        END IF
        CALL CUTEST_probname( status, pname )
        WRITE( summary_unit, "( A10 )" ) pname
      END IF

      IF ( iter_summary_unit > 0 .and. control%output_progress_vectors ) THEN
        INQUIRE( FILE = iter_summary_file, EXIST = filexx )
        IF ( filexx ) THEN
           OPEN( iter_summary_unit, FILE=iter_summary_file, FORM='FORMATTED',  &
               STATUS = 'OLD', IOSTAT = iores , position="append")
        ELSE
           OPEN( iter_summary_unit, FILE=iter_summary_file, FORM='FORMATTED',  &
                STATUS = 'NEW', IOSTAT = iores )
        END IF
        IF ( iores /= 0 ) THEN 
          write( out, "( ' IOSTAT = ', I0, ' when opening file ', A,           &
        &  '. Stopping ' )" ) iores, iter_summary_file
          STOP
        END IF
        CALL CUTEST_probname( status, pname )
        WRITE( iter_summary_unit, "( A10 )" ) pname
      END IF
      
      sol_unit = 17
      open( sol_unit, FILE="last.solution", FORM='FORMATTED',   &
           STATUS= 'REPLACE', IOSTAT = iores)!, position='REWIND' )

      !  set up structures for eval_hp
      call cutest_cdimchp(status, params%lchp)
      IF ( status /= 0) GO TO 910
      allocate(params%chp_ptr(m+1))
      allocate(params%chp_ind(params%lchp))
      allocate(params%chp_val(params%lchp))
      params%GotH = .false.
      
      
      write(*,*) 'calling the minimizer...'

!  call the minimizer
      write(*,*) 'sending to nlls_solve'
      if (supply_eval_hp == 1) then 
         CALL NLLS_SOLVE( n, m, X, eval_F, eval_J, eval_HF,                         &
              params, control, inform, eval_HP=eval_HP )
      else
         CALL NLLS_SOLVE( n, m, X, eval_F, eval_J, eval_HF,                         &
              params, control, inform)
      end if
!           params, control, inform)

      WRITE( out , "( A, I0, A, I0)") 'status = ', inform%status,              &
          '       iter = ', inform%iter
      IF ( status /= 0 ) GO TO 910

!  output report

      CALL CUTEST_creport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910

      ALLOCATE( F( m ), VNAMES( n ), CNAMES( m ), STAT = status )
      CALL CUTEST_cnames( status, n, m, pname, VNAMES, CNAMES )
      CALL eval_F( status, n, m, X, F, params)

      WRITE( out, 2110 ) ( i, VNAMES( i ), X( i ), i = 1, n )
!     WRITE( out, 2120 ) ( i, CNAMES( i ), F( i ), i = 1, m )
      WRITE( out, 2000 ) pname, n, m, inform%obj, INT( CALLS( 5 ) ),           &
        INT( CALLS( 6 ) ), INT( CALLS( 7 ) ), CPU( 1 ), CPU( 2 )

!  write summary if required

      localiter = inform%iter
      fnevals = inform%f_eval!int( calls(5) ) 
      jacevals = inform%g_eval!int( calls(6) ) 
      hessevals = inform%h_eval!int( calls(7) ) 
      inner_iter = inform%inner_iter
      solve_time = CPU(2)
      
      if ( inform%status .ne. 0 ) then 
         localiter = -localiter
         fnevals = -fnevals
         jacevals = -jacevals
         hessevals = -hessevals
      end if

      IF ( summary_unit > 0 ) THEN
        BACKSPACE( summary_unit )
        WRITE( summary_unit, "( A10, 7I8, I8, 4ES25.15E3 )" ) & !, ES23.15E3, ES23.15E3, ES23.15E3 )" ) &
          pname, n, m, inform%status,localiter,                                &
          fnevals, jacevals, hessevals,                                        &
          inner_iter,                                                          &
          inform%obj, inform%norm_g, inform%scaled_g, abs(solve_time)
        CLOSE(  summary_unit )
      END IF
      
      IF ( iter_summary_unit > 0 .and. control%output_progress_vectors ) THEN
        BACKSPACE( iter_summary_unit )
        do i = 1,inform%iter + 1
           WRITE( iter_summary_unit, "( ES23.15E3,A, ES23.15E3 )" )              &
                inform%resvec(i),'  ', inform%gradvec(i)
        end do
        CLOSE(  iter_summary_unit )
      END IF
      do i = 1, n
         WRITE( sol_unit, "( ES23.15E3 )") X(i)
      end do
!  clean-up data structures
      DEALLOCATE( X, F, VNAMES, CNAMES, STAT = status )
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

2000 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //,                    &
          ' Package used            :  RAL_NLLS ',  /,                         &
          ' Problem                 :  ', A10,    /,                           &
          ' # variables             =  ', I0, /,                               &
          ' # residuals             =  ', I0, /,                               &
          ' Final f                 =', ES15.7 /,                              &
          ' # residual evaluations  =  ', I0, /,                               &
          ' # Jacobian evaluations  =  ', I0, /,                               &
          ' # Hessian evaluations   =  ', I0, /,                               &
          ' Set up time             =  ', 0P, F0.2, ' seconds' /,              &
          ' Solve time              =  ', 0P, F0.2, ' seconds' //,             &
           66('*') / )
2110 FORMAT( /, ' The variables:', /, &
          '     i name          value',  /, ( I6, 1X, A10, 1P, D12.4 ) )
!2120 FORMAT( /, ' The constraints:', /, '     i name          value',         &
!         /, ( I6, 1X, A10, 1P, D12.4 ) )

!  End of RAL_NLLS_main

    contains

      SUBROUTINE eval_F( status, n, m, X, F, params )
      use :: ral_nlls_double, only : params_base_type
      integer, intent(out) :: status
      integer, intent(in) :: n,m
      double precision, dimension(*), intent(in)  :: x
      double precision, dimension(*), intent(out) :: f
      class(params_base_type), intent(inout) :: params
      double precision :: obj
!  evaluate the residuals F

      CALL CUTEST_cfn( status, n, m, X, obj, F )
      END SUBROUTINE eval_F

      SUBROUTINE eval_J( status, n, m, X, J, params)
      USE ISO_C_BINDING
      use :: ral_nlls_double, only : params_base_type
      
      INTEGER ( c_int ), INTENT( OUT ) :: status
      INTEGER ( c_int ), INTENT( IN ) :: n, m
      REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
      REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: J
      class( params_base_type ), intent(inout) :: params
      REAL ( c_double ), DIMENSION( n ) :: G
      REAL ( c_double ), DIMENSION( m ) :: Y
      REAL ( c_double ), DIMENSION( m , n ) :: Jmatrix

!  evaluate the residual Jacobian J

      CALL CUTEST_cgr( status, n, m, X, Y, .FALSE., G, .FALSE., m, n, Jmatrix ) 
      ! convert the Jacobian to a vector....
      J(1:m*n) = reshape(Jmatrix, (/n*m/) )
      RETURN
      END SUBROUTINE eval_J

      SUBROUTINE eval_HF( status, n, m, X, F, H, params )
        USE ISO_C_BINDING
        use :: ral_nlls_double, only : params_base_type
        
        INTEGER ( c_int ), INTENT( OUT ) :: status
        INTEGER ( c_int ), INTENT( IN ) :: n, m
        REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: X
        REAL ( c_double ), DIMENSION( * ), INTENT( IN ) :: F
        REAL ( c_double ), DIMENSION( * ), INTENT( OUT ) :: H
        class( params_base_type ), intent(inout) :: params
        
        real ( c_double ), dimension(n,n) :: Hmatrix
        !  evaluate the product H = sum F_i Hessian F_i

        CALL CUTEST_cdhc( status, n, m, X, F, n, Hmatrix )
        
        H(1:n*n) = reshape(Hmatrix, (/n*n/) )
        RETURN
      END SUBROUTINE eval_HF

      subroutine eval_HP(status, n, m, x, y, hp, params)
        use iso_c_binding
        use :: ral_nlls_double, only : params_base_type
        
        integer ( c_int ), intent( out ) :: status
        integer ( c_int ), intent( in ) :: n, m
        real ( c_double ), dimension(*), intent( in ) :: x, y
        real ( c_double ), dimension(*), intent( out ) :: hp
        class( params_base_type ), intent(inout) :: params

        integer :: i, j

        select type(params)
        type is(user_type)           
           call cutest_cchprods(status, n, m, params%GotH, x, y, params%lchp, &
                params%chp_val(1:params%lchp), params%chp_ind(1:params%lchp), params%chp_ptr(1:m+1))
           hp(1:n*m) = 0.0
           do i =  1,m ! loop over the columns
              do j = params%chp_ptr(i),params%chp_ptr(i+1)-1 ! and the rows...
                 hp((i-1)*n + params%chp_ind(j)) = params%chp_val(j)
              end do
           end do
           params%GotH = .true.
        end select
        
      end subroutine eval_HP

    END PROGRAM RAL_NLLS_main





