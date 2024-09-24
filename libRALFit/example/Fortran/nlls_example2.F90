! Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
! All rights reserved.
! Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
! examples/Fortran/nlls_example2.f90
!
! Attempts to fit the model y_i = x_1 e^(x_2 t_i)
! For parameters x_1 and x_2, and input data (t_i, y_i)
module fndef_example2

#if SINGLE_PRECISION
   use ral_nlls_single
#else
   use ral_nlls_double
#endif

   implicit none

   type, extends(params_base_type) :: params_type
      real(wp), dimension(:), allocatable :: t ! The m data points t_i
      real(wp), dimension(:), allocatable :: y ! The m data points y_i
   end type

contains
   ! Calculate r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
   subroutine eval_r(status, n, m, x, r, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: r
      class(params_base_type), intent(inout) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(n)
      select type(params)
      type is(params_type)
         r(1:m) = x1 * exp(x2*params%t(:)) - params%y(:)
      end select

      status = 0 ! Success
   end subroutine eval_r
   ! Calculate:
   ! J_i1 = e^(x_2 * t_i)
   ! J_i2 = t_i x_1 e^(x_2 * t_i)
   subroutine eval_J(status, n, m, x, J, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: J
      class(params_base_type), intent(inout) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(n)
      select type(params)
      type is(params_type)
         J(  1:  m) = exp(x2*params%t(1:m))                     ! J_i1
         J(m+1:2*m) = params%t(1:m) * x1 * exp(x2*params%t(1:m))! J_i2
      end select

      status = 0 ! Success
   end subroutine eval_J
   ! Calculate
   ! HF = sum_i r_i H_i
   ! Where H_i = [ 0                t_i e^(x_2 t_i)    ]
   !             [ t_i e^(x_2 t_i)  t_i^2 x_1 e^(x_2 t_i)  ]
   subroutine eval_HF(status, n, m, x, r, HF, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(in) :: r
      real(wp), dimension(*), intent(out) :: HF
      class(params_base_type), intent(inout) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(2)
      select type(params)
      type is(params_type)
         HF(    1) = sum(r(1:m) * 0)                                      ! H_11
         HF(    2) = sum(r(1:m) * params%t(1:m) * exp(x2*params%t(1:m)))  ! H_21
         HF(1*n+1) = HF(2)                                                ! H_12
         HF(1*n+2) = sum(r(1:m) * (params%t(1:m)**2) * x1 * exp(x2*params%t(1:m)))! H_22
      end select

      status = 0 ! Success
    end subroutine eval_HF

    subroutine eval_HP(status, n, m, x, y, HP, params)
      implicit none
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(in) :: y
      real(wp), dimension(*), intent(out) :: HP
      class(params_base_type), intent(inout) :: params

      real(wp) :: x1, x2
      integer :: i

      x1 = x(1)
      x2 = x(2)
      select type(params)
      type is(params_type)
         do i = 1, m
            HP(n*(i-1) + 1) = 0.0*y(1) + params%t(i)*exp(x2*params%t(i))*y(2)
            HP(n*(i-1) + 2) = params%t(i) * exp( x2*params%t(i) ) * y(1) + &
                              (params%t(i)**2) * x1 * exp(x2*params%t(i))*y(2)
         end do
      end select

      status = 0 ! Success
    end subroutine eval_HP

end module fndef_example2

program nlls_example2
   use fndef_example2
   implicit none

   type(nlls_options) :: options
   type(nlls_inform) :: inform

   integer :: m,n
   real(wp), allocatable :: x(:)
   type(params_type) :: params
   integer :: inner_method
   ! Data to be fitted
   m = 5
   allocate(params%t(m), params%y(m))
   params%t(:) = (/ 1.0, 2.0, 4.0,  5.0,  8.0 /)
   params%y(:) = (/ 3.0, 4.0, 6.0, 11.0, 20.0 /)

   ! Call fitting routine
   n = 2
   allocate(x(n))

   options%model = 4
   options%type_of_method = 2
   options%nlls_method = 4
   options%exact_second_derivatives = .true.
   options%maxit = 50
   options%print_level = 4

   do inner_method = 1,3

      options%inner_method = inner_method

      x = (/ 2.5, 0.25 /) ! Initial guess
      call nlls_solve(n, m, x, eval_r, eval_J, eval_HF, params, options, inform)
      if(inform%status.ne.0) then
         print *, "ral_nlls() returned with error flag ", inform%status
         goto 100
      endif

      ! Print result
      print *, "Found a local optimum at x = ", x
      print *, "Took ", inform%iter, " iterations"
      print *, "     ", inform%f_eval, " function evaluations"
      print *, "     ", inform%g_eval, " gradient evaluations"
      print *, "     ", inform%h_eval, " hessian evaluations"
!!$
!!$      print *, "===== passing eval_HP ====="
!!$
!!$      x = (/ 2.5, 0.25 /) ! Initial guess
!!$      call nlls_solve(n, m, x, eval_r, eval_J, eval_HF, params, options, inform,eval_HP=eval_HP)
!!$      if(inform%status.ne.0) then
!!$         print *, "ral_nlls() returned with error flag ", inform%status
!!$         goto 100
!!$      endif
!!$
!!$      ! Print result
!!$      print *, "Found a local optimum at x = ", x
!!$      print *, "Took ", inform%iter, " iterations"
!!$      print *, "     ", inform%f_eval, " function evaluations"
!!$      print *, "     ", inform%g_eval, " gradient evaluations"
!!$      print *, "     ", inform%h_eval, " hessian evaluations"

   end do

100 Continue

   if (allocated(x)) deallocate(x)
   if (allocated(params%t)) deallocate(params%t)
   if (allocated(params%y)) deallocate(params%y)

end program nlls_example2
