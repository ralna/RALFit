! Interface to DCT implementations.
#include "../preprocessor.FPP"

module MODULE_PREC(dct_module)
  
  implicit none

#ifdef SINGLE_PRECISION
  integer, parameter :: wp = selected_real_kind(6) ! c_float
#else
  integer, parameter :: wp = selected_real_kind(15) ! c_double
#endif

  ! these parameters are defined in FFTW, but we define them here to avoid including the FFTW header
  integer :: FFTW_FORWARD 
  parameter (FFTW_FORWARD = -1)
  integer :: FFTW_ESTIMATE
  parameter (FFTW_ESTIMATE = 64)

  interface dct
    module procedure dct_fftz
  end interface dct

  interface dct1d
    module procedure dct_fftz_1d
  end interface dct1d

  !interface fftw_plan_dft_1d
  !  module subroutine fftw_plan_dft_1d(plan, n, in, out, sign, flags) bind(C)
  !    use iso_c_binding, only: c_int, c_float, c_double
  !    integer(c_int) :: plan
  !    integer(c_int), value :: n
  !    complex(wp) :: in(*)
  !    complex(wp) :: out(*)
  !    integer(c_int), value :: sign
  !    integer(c_int), value :: flags
  !  end subroutine fftw_plan_dft_1d
  !end interface fftw_plan_dft_1d
  
  contains

    subroutine dct_fftz(A, n, m, work)
      ! Compute the 1D DCT across the columns of A using the FFTZ library.
      ! `work` is expected to be a complex array of size m.
      real(wp), intent(inout), contiguous :: A(:,:)
      complex(wp), intent(inout) :: work(:)
      integer, intent(in) :: n, m 

      integer :: i, plan

      plan = 0

      do i = 1, n
        call dct_fftz_1d(A(:, i), m, work)
      end do 

    end subroutine dct_fftz

    subroutine dct_fftz_1d(x, m, work)
      ! Compute the 1D DCT of x using the FFTZ library.
      ! This uses an FFT and a phase factor to compute the DCT coefficients,
      ! because it's easier to get FFT codes than DCT ones.
      ! `work` is expected to be a complex array of size m, 
      ! and `plan` is an integer that can be reused across calls 
      ! to avoid the overhead of creating a new plan each time.
      real(wp), intent(inout), contiguous :: x(:)
      complex(wp), intent(inout) :: work(:)
      complex(wp) :: phase_factor
      integer, intent(in) :: m

      integer :: i, plan

      ! phase factor for converting FFT output into DCT coefficients
      phase_factor = complex(0.0_wp, -3.14159_wp/(2.0_wp*m))

      ! The DCT can be computed using an FFT by creating a new array where the first m/2 elements are the even-indexed
      ! elements of x and the next m/2 elements are the odd-indexed elements of x (in reverse). 
      ! Then we can compute the DCT using an FFT on this new array.
      do i = 1, m/2
        work(i) = complex(x(2*i - 1), 0.0_wp)
        work(m/2 + i) = complex(x(2*(m/2 - i + 1)), 0.0_wp)
      end do

      call dfftw_plan_dft_1d(plan, m, work, work, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan, work, work)

      do i = 1, m
        x(i) = real(work(i) * exp(phase_factor * (i-1)), wp)
      end do
      
      call dfftw_destroy_plan(plan)

    end subroutine dct_fftz_1d

end module MODULE_PREC(dct_module)
