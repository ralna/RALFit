! THIS VERSION: RAL_NLLS 1.0 - 22/12/2015 AT 14:15 GMT.

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*     SYMBOLS  M O D U L E    *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   extracted from GALAHAD package SYMBOLS, December 22nd, 2015

#include "preprocessor.FPP"

   MODULE MODULE_PREC(RAL_NLLS_SYMBOLS)

!  This module provides the list of all symbolic names that are common to
!  the RAL_NLLS modules. It is just intended as a dictionary for use in other
!  modules.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      IMPLICIT NONE
      PUBLIC :: SYMBOLS_status

!  exit conditions (0 to -99; others will be package specific)

      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_ok                      = 0
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_allocate          = - 1
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_deallocate        = - 2
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_restrictions      = - 3
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_bad_bounds        = - 4
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_primal_infeasible = - 5
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_dual_infeasible   = - 6
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unbounded         = - 7
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_no_center         = - 8
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_analysis          = - 9
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_factorization     = - 10
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_solve             = - 11
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_analysis      = - 12
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_factorization = - 13
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_uls_solve         = - 14
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_preconditioner    = - 15
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ill_conditioned   = - 16
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_tiny_step         = - 17
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_max_iterations    = - 18
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_time_limit        = - 19
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cpu_limit         =         &
                                    RAL_NLLS_error_time_limit
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_inertia           = - 20
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_file              = - 21
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_io                = - 22
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_upper_entry       = - 23
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_sort              = - 24
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_input_status      = - 25
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unknown_solver    = - 26
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_not_yet_implemented     = - 27
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qp_solve          = - 28
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_unavailable_option      = - 29
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_warning_on_boundary     = - 30
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_call_order        = - 31
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_integer_ws        = - 32
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_real_ws           = - 33
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_pardiso           = - 34
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_wsmp              = - 35
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc64              = - 36
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc77              = - 37
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_lapack            = - 38
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_permutation       = - 39
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_alter_diagonal    = - 40
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_access_pivots     = - 41
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_access_pert       = - 42
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_direct_access     = - 43
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_f_min             = - 44
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_unknown_precond   = - 45
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_schur_complement  = - 46
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_technical         = - 50
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_reformat          = - 52
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ah_unordered      = - 53
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_y_unallocated     = - 54
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_z_unallocated     = - 55
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_scale             = - 61
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_presolve          = - 62
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpa               = - 63
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpb               = - 64
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_qpc               = - 65
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cqp               = - 66
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_dqp               = - 67
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc61              = - 69
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mc68              = - 70
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_metis             = - 71
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_spral             = - 72
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_warning_repeated_entry  = - 73
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_rif               = - 74
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ls28              = - 75
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ls29              = - 76
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_cutest            = - 77
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_evaluation        = - 78
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_optional          = - 79
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_mi35              = - 80
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_spqr              = - 81
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_alive             = - 82
      INTEGER, PUBLIC, PARAMETER :: RAL_NLLS_error_ccqp              = - 83

   CONTAINS

!-*-  R A L _ N L L S -  S Y M B O L S _ S T A T U S   S U B R O U T I N E  -*-

     SUBROUTINE SYMBOLS_status( status, out, prefix, routine )

!  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!  Print details of return codes

!  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: status, out
     CHARACTER ( LEN = * ) :: prefix
     CHARACTER ( LEN = * ) :: routine

     SELECT CASE ( status )
     CASE( RAL_NLLS_ok )

     CASE( RAL_NLLS_error_allocate )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   allocation error' )" )                                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_deallocate )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   deallocation error' )" )                                     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_restrictions )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
     &        '   one or more restrictions violated' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_bad_bounds )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem bounds are infeasible' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_primal_infeasible )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem appears to be (locally) infeasible' )" )         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_dual_infeasible )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the problem appears to be (locally) dual infeasible' )" )    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unbounded )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
     &        '   the problem appears to be unbounded from below' )" )         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_no_center )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the feasible region has no analytic center' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_analysis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the analysis phase failed' )" )                              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_factorization )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the factorization failed' )" )                               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the solve phase failed' )" )                                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_analysis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &        '   the unsymmetric analysis phase failed' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_factorization )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the unsymmetric factorization failed' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_uls_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the unsymmetric solve phase failed' )" )                     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_preconditioner )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the provided preconditioner is flawed' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ill_conditioned )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A, '   the',   &
      &       ' problem is too ill conditioned to make further progress' )" )  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_tiny_step )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A, '   the',   &
      &       ' computed step is too small to make further progress' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_max_iterations )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the iteration limit has been exceeded' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_time_limit )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the time limit has been exceeded' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_inertia )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the preconditioner has the wrong inertia' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_file )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a file-handling error occurred' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_io )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an input/output error occurred' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_upper_entry )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   there is a matrix entry in the upper triangle' )" )          &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_sort )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occurred when sorting' )" )                         &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_input_status )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   error with input status' )" )                                &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unknown_solver )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested solver is not known' )" )                      &
         prefix, routine, prefix
     CASE ( RAL_NLLS_not_yet_implemented )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested option has not yet been implemented' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qp_solve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   QP solver failed: check the QP solver status' )" )           &
         prefix, routine, prefix
     CASE( RAL_NLLS_unavailable_option )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested option is unavailable' )" )                    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_call_order )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the routine has been called out of order' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_integer_ws )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   insufficient integer workspace for factorization' )" )       &
         prefix, routine, prefix, status
     CASE( RAL_NLLS_error_real_ws )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   insufficient real workspace for factorization' )" )          &
         prefix, routine, prefix, status
     CASE( RAL_NLLS_error_pardiso )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   PARDISO failure: check its return status' )" )               &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_wsmp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   WSMP failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc64 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC64 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc77 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC77 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_lapack )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LAPACK failure: check its return status' )" )                &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_permutation )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   input order is not a permutation' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_alter_diagonal )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to alter diagonals with this solver' )" )   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_access_pivots )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to access the pivots with this solver' )" ) &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_access_pert )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   it is impossible to access the perturbations made with',     &
      &       ' this solver' )" ) prefix, routine, prefix
     CASE( RAL_NLLS_error_direct_access )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a direct-access file error occurred' )" )                    &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_f_min )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the objective value is too small' )" )                       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_unknown_precond )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested preconditioner is not known' )" )              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_schur_complement )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the requested Schur complement is too large' )" )            &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_technical )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   a technical error has occurred within the linear solver' )" )&
         prefix, routine, prefix
     CASE( RAL_NLLS_error_reformat )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the matrix storage format has been changed' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ah_unordered )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   the matrix storage format should have been changed' )" )     &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_y_unallocated )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   Y, Y_l or Y_u has not been allocated' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_z_unallocated )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   Z, Z_l or Z_u has not been allocated' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_scale )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when scaling the problem' )" )              &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_presolve )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when pre/postsolving the problem' )" )      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpa )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPA' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpb )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPB' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_qpc )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling QPC' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_cqp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling CQP' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_dqp )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   an error occured when calling DQP' )" )                      &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc61 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC61 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mc68 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MC68 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_metis )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   METIS failure: check its return status' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_spral )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   SPRAL failure: check its return status' )" )                 &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ls28 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LS28 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_ls29 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   LS29 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_mi35 )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   MI35 failure: check its return status' )" )                  &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_rif )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   RIF failure: check its return status' )" )                   &
         prefix, routine, prefix
     CASE( RAL_NLLS_warning_repeated_entry )
       WRITE( out, "( /, A,  ' Warning return from ', A, ' -', /, A,           &
      &       '   the input matrix contains repeated entries' )" )             &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_cutest )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   CUTEst evaluation error: check its return status' )" )       &
         prefix, routine, prefix
     CASE( RAL_NLLS_error_evaluation )
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   function evaluation error: check its return status' )" )     &
         prefix, routine, prefix
     CASE DEFAULT
       WRITE( out, "( /, A,  ' Error return from ', A, ' -', /, A,             &
      &       '   status = ', I0 )" ) prefix, routine, prefix, status
     END SELECT

     RETURN

!  End of subroutine SYMBOLS_status

     END SUBROUTINE SYMBOLS_status

   END MODULE MODULE_PREC(RAL_NLLS_SYMBOLS)
