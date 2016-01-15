! examples/Fortran/nlls_example.f90
! 
! Attempts to fit the model y_i = x_1 e ^(x_2 t_i)
! For parameters x_1 and x_2, and input data (t_i, y_i)
module fndef_example
   use nlls_module, only : params_base_type
   implicit none

   integer, parameter :: wp = kind(0d0)

   type, extends(params_base_type) :: params_type
      real(wp), dimension(:), allocatable :: t ! The m data points t_i
      real(wp), dimension(:), allocatable :: y ! The m data points y_i
   end type

contains
   ! Calculate f_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
   subroutine eval_F(status, n, m, x, f, params)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: f
      class(params_base_type), intent(in) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(2)
      select type(params)
      type is(params_type)
         f(1:m) = x1 * exp(x2*params%t(:)) - params%y(:)
      end select

      status = 0 ! Success
   end subroutine eval_F
   ! Calculate:
   ! J_i1 = e^(x_2 * t_i)
   ! J_i2 = t_i x_1 e^(x_2 * t_i)
   subroutine eval_J(status, n, m, x, J, params)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(out) :: J
      class(params_base_type), intent(in) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(2)
      select type(params)
      type is(params_type)
         J(  1:  m) = exp(x2*params%t(:))                      ! J_i1
         J(m+1:2*m) = params%t(:) * x1 * exp(x2*params%t(:))   ! J_i2
      end select

      status = 0 ! Success
   end subroutine eval_J
   ! Calculate
   ! f_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
   ! HF = sum_i f_i H_i
   ! Where H_i = [ 1                t_i e^(x_2 t_i)    ]
   !             [ t_i e^(x_2 t_i)  t_i^2 e^(x_2 t_i)  ]
   subroutine eval_HF(status, n, m, x, f, HF, params)
      integer, intent(out) :: status
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(wp), dimension(*), intent(in) :: x
      real(wp), dimension(*), intent(in) :: f
      real(wp), dimension(*), intent(out) :: HF
      class(params_base_type), intent(in) :: params

      real(wp) :: x1, x2

      x1 = x(1)
      x2 = x(2)
      select type(params)
      type is(params_type)
         HF(1) = sum(f(1:m) * 1)                                     ! H_11
         HF(2) = sum(f(1:m) * params%t(:) * exp(x2*params%t(:)))     ! H_21
         HF(3) = HF(2)                                               ! H_12
         HF(4) = sum(f(1:m) * params%t(:)**2 * exp(x2*params%t(:)))  ! H_22
      end select

      status = 0 ! Success
   end subroutine eval_HF
end module fndef_example

program nlls_example
   use nlls_module
   use fndef_example
   implicit none

   type(nlls_control_type) :: control
   type(nlls_inform_type) :: info

   integer :: m
   real(wp), dimension(2) :: x
   type(params_type) :: params

   ! Data to be fitted
   m = 5
   allocate(params%t(5), params%y(5))
   params%t(:) = (/ 1.0, 2.0, 4.0,  5.0,  8.0 /)
   params%y(:) = (/ 3.0, 4.0, 6.0, 11.0, 20.0 /)

   ! Call fitting routine
   control%nlls_method = 4 ! FIXME: remove
   x(:) = (/ 2.5, 0.25 /) ! Initial guess
   call ral_nlls(2, m, x, eval_F, eval_J, eval_HF, params, info, control)
   if(info%status.ne.0) then
      print *, "ral_nlls() returned with error flag ", info%status
      stop
   endif

   ! Print result
   print *, "Found a local optimum at x = ", x
   print *, "Residuals = ", info%resvec(:)
   print *, "Took ", info%iter, " iterations"
   print *, "     ", info%f_eval, " function evaluations"
   print *, "     ", info%g_eval, " gradient evaluations"
   print *, "     ", info%h_eval, " hessian evaluations"
end program nlls_example
