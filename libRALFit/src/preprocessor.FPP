! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.

#define CONCAT_(a) a
#define CONCAT(a, b) CONCAT_(a)CONCAT_(b)

#if SINGLE_PRECISION

! BLAS nomenclature
#define PREC(fn) CONCAT(s,fn)
! Module precion nomenclature
#define MODULE_PREC(fn) CONCAT(fn,_single)
! CIFACE nomenclature
#define IFACE_PREC(fn) CONCAT(fn,_s)

#else

    ! BLAS nomenclature
#define PREC(fn) CONCAT(d,fn)
! Module precion nomenclature
#define MODULE_PREC(fn) CONCAT(fn,_double)
! CIFACE nomenclature
#define IFACE_PREC(fn) CONCAT(fn,_d)

#endif