! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_types :: Globally defines the types used
!
! Working precision (wp) for reals(Kind=wp)
! exports for C binding:
! integer type (ILP64) and aliases (ral_c_int => default integer size (4 or 8)
! real type float/double aliases to wp: ral_c_real => wp

#include "preprocessor.FPP"

module MODULE_PREC(ral_nlls_types)

   use iso_c_binding, only : c_int, c_int64_t, c_double, c_float

   implicit none

   Public

   Integer, Parameter :: np = selected_real_kind(15) ! c_double
   Integer, Parameter :: lp = selected_real_kind(6) ! c_float

#ifdef SINGLE_PRECISION
   Integer, Parameter :: wp = lp
   Integer, Parameter :: ral_c_real = c_float
#else
   Integer, Parameter :: wp = np
   Integer, Parameter :: ral_c_real = c_double
#endif

#ifdef BUILD_ILP64
   ! To build ilp64 also the default integer kind must be set to 8 bytes
   Integer, Parameter :: ral_c_int = c_int64_t
#else
   Integer, Parameter :: ral_c_int = c_int
#endif

  type :: params_base_type
     ! deliberately empty
  end type params_base_type

! Abstract interfaces
  abstract interface
     subroutine eval_f_type(status, n, m, x, f, params)
       import :: wp,params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       Real(Kind=wp), dimension(*), intent(in)  :: x
       Real(Kind=wp), dimension(*), intent(out) :: f
       class(params_base_type), intent(inout) :: params
     end subroutine eval_f_type
  end interface

  abstract interface
     subroutine eval_j_type(status, n, m, x, J, params)
       import :: wp,params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       Real(Kind=wp), dimension(*), intent(in)  :: x
       Real(Kind=wp), dimension(*), intent(out) :: J
       class(params_base_type), intent(inout) :: params
     end subroutine eval_j_type
  end interface
  abstract interface
     subroutine eval_hf_type(status, n, m, x, f, h, params)
       import :: wp,params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       Real(Kind=wp), dimension(*), intent(in)  :: x
       Real(Kind=wp), dimension(*), intent(in)  :: f
       Real(Kind=wp), dimension(*), intent(out) :: h
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hf_type
  end interface

  abstract interface
     subroutine eval_hp_type(status, n, m, x, y, hp, params)
       import :: wp,params_base_type
       implicit none
       integer, intent(out) :: status
       integer, intent(in) :: n,m
       Real(Kind=wp), dimension(*), intent(in)  :: x
       Real(Kind=wp), dimension(*), intent(in)  :: y
       Real(Kind=wp), dimension(*), intent(out) :: hp
       class(params_base_type), intent(inout) :: params
     end subroutine eval_hp_type
  end interface


! Constants related to types
#ifdef SINGLE_PRECISION
   ! max display largest (2 digit exponent)
   Real(Kind=wp), Parameter :: prn_pinf = +9.99e+37_wp
   Real(Kind=wp), Parameter :: prn_minf = -9.99e+37_wp
   Real(Kind=wp), Parameter :: prn_small = 1.00e-37_wp

   ! constants related to testing
   Real(Kind=wp), Parameter :: test_big = +1.00e+37_wp
   Real(Kind=wp), Parameter :: test_very_big = +2.00e+37_wp
   Real(Kind=wp), Parameter :: test_very_small = 1.00e-37_wp
#else
   ! max display largest (2 digit exponent)
   Real(Kind=wp), Parameter :: prn_pinf = +9.99e+99_wp
   Real(Kind=wp), Parameter :: prn_minf = -9.99e+99_wp
   Real(Kind=wp), Parameter :: prn_small = 1.00e-99_wp

   ! constants related to testing
   Real(Kind=wp), Parameter :: test_big = +1.00e+100_wp
   Real(Kind=wp), Parameter :: test_very_big = +1.00e+105_wp
   Real(Kind=wp), Parameter :: test_very_small = 1.00e-100_wp
#endif

end module MODULE_PREC(ral_nlls_types)

