! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! ral_nlls_types :: Globally defines the types used
!
! Working precision (wp) for reals(Kind=wp)
! exports for C binding:
! integer type (ILP64) and aliases (ral_c_int => ip
! real type float/double aliases to wp: ral_c_real => wp

module ral_nlls_types

   use iso_c_binding, only : c_int, c_int64_t, c_double, c_float

   implicit none

   Integer, Parameter, Public :: np = selected_real_kind(15) ! c_double
   Integer, Parameter, Public :: lp = selected_real_kind(6) ! c_float


#ifdef SINGLE_PRECISION
   Integer, Parameter, Public :: wp = lp
   Integer, Parameter, Public :: ral_c_real = c_float
#else
   Integer, Parameter, Public :: wp = np
   Integer, Parameter, Public :: ral_c_real = c_double
#endif

#ifdef BUILD_ILP64
   ! To build ilp64 also the default integer kind must be set to 8 bytes
   Integer, Parameter, Public :: ral_c_int = c_int64_t
#else
   Integer, Parameter, Public :: ral_c_int = c_int
#endif

end module ral_nlls_types

