! Interface to DCT implementations.
#include "../preprocessor.FPP"

module MODULE_PREC(dct_module)
  implicit none

  ! these parameters are defined in FFTW, but we define them here to avoid including the FFTW header
  integer(8) :: FFTW_REDFT10
  parameter (FFTW_REDFT10 = 4)
  integer(8) :: FFTW_ESTIMATE
  parameter (FFTW_ESTIMATE = 64)
  integer(4),  parameter :: wp = kind( 0.0d+0 )

  interface dct
    module procedure dct_fftz
  end interface dct
  
  contains

    subroutine dct_fftz(A, n, m)
      ! Compute the 1D DCT across the columns of A using the FFTZ library.
      real(8), intent(inout), contiguous :: A(:,:)
      integer, intent(in) :: n, m 

      integer :: i, plan

      call dfftw_plan_r2r_1d(plan, m, A(:, 1), A(:, 1), FFTW_REDFT10, FFTW_ESTIMATE)
      do i = 1, n
        call dfftw_execute_r2r(plan, A(:, i), A(:, i))
      end do 
      call dfftw_destroy_plan(plan)

    end subroutine dct_fftz

end module MODULE_PREC(dct_module)
