module example_module

  use :: ral_nlls_double, only : params_base_type
  implicit none 

  type, extends( params_base_type ) :: user_type
     double precision, allocatable :: x_values(:)
     double precision, allocatable :: y_values(:)
  end type user_type

contains
  
  
SUBROUTINE eval_F( status, n, m, X, f, params)

!  -------------------------------------------------------------------
!  eval_F, a subroutine for evaluating the function f at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER ( int32 ), INTENT( OUT ) :: status
       INTEGER ( int32 ), INTENT( IN ) :: n, m 
       REAL ( real64 ), DIMENSION( * ),INTENT( OUT ) :: f
       REAL ( real64 ), DIMENSION( * ),INTENT( IN )  :: X
       class( params_base_type ), intent(in) :: params
! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i

! then, let's work this into the format we need
! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          do i = 1,m
             f(i) = params%y_values(i) - exp( X(1) * params%x_values(i) + X(2) )
          end do
       end select

       status = 0
       
!!$! let's use Powell's function for now....
!!$       f(1) = X(1) + 10.0 * X(2)
!!$       f(2) = sqrt(5.0) * (X(3) - X(4))
!!$       f(3) = ( X(2) - 2.0 * X(3) )**2
!!$       f(4) = sqrt(10.0) * ( X(1) - X(4) )**2
       
! end of subroutine eval_F
       
     END SUBROUTINE eval_F


     SUBROUTINE eval_J( status, n, m, X, J, params)

!  -------------------------------------------------------------------
!  eval_J, a subroutine for evaluating the Jacobian J at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER ( int32 ), INTENT( OUT ) :: status
       INTEGER ( int32 ), INTENT( IN ) :: n, m 
       REAL ( real64 ), DIMENSION( * ),INTENT( OUT ) :: J
       REAL ( real64 ), DIMENSION( * ),INTENT( IN ) :: X
       class( params_base_type ), intent(in) :: params

! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i

       ! let's work this into the format we need
       ! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          do i = 1,m
             J(i) =  - params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) )
             J(m + i) = - exp( X(1) * params%x_values(i) + X(2) )
          end do
       end select
       
       status = 0

! end of subroutine eval_J

!!$       ! initialize to zeros...
!!$       J(1:4,1:4) = 0.0
!!$       
!!$       ! enter non-zeros values
!!$       J(1,1) = 1.0
!!$       J(1,2) = 10.0
!!$       J(2,3) = sqrt(5.0)
!!$       J(2,4) = -sqrt(5.0)
!!$       J(3,2) = 2.0 * (X(2) - 2.0 * X(3))
!!$       J(3,3) = -4.0 * (X(2) - 2.0 * X(3)) 
!!$       J(4,1) = sqrt(10.0) * 2.0 * (X(1) - X(4))
!!$       J(4,4) = - sqrt(10.0) * 2.0 * (X(1) - X(4))

     END SUBROUTINE eval_J


SUBROUTINE eval_H( status, n, m, X, f, h, params)

!  -------------------------------------------------------------------
!  eval_F, a subroutine for evaluating the function f at a point X
!  -------------------------------------------------------------------

       USE ISO_FORTRAN_ENV

       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER ( int32 ), INTENT( OUT ) :: status
       INTEGER ( int32 ), INTENT( IN ) :: n, m 
       REAL ( real64 ), DIMENSION( * ),INTENT( IN )  :: f
       REAL ( real64 ), DIMENSION( * ),INTENT( OUT ) :: h
       REAL ( real64 ), DIMENSION( * ),INTENT( IN )  :: X
       class( params_base_type ), intent(in) :: params
! Let's switch to an actual fitting example...
! min 0.5 || f(m,c)||**2, where
! f_i(m,c) = y_i - exp( m * x_i + c )

       integer :: i

! then, let's work this into the format we need
! X(1) = m, X(2) = c
       select type(params)
       type is(user_type)
          ! evaluate 
          ! HF = \sum_{i=1}^m F_i H_i
          h(1:4) = 0.0
          do i = 1, m
             h(1) = &
                  h(1) + f(i)* ( & 
                  - (params%x_values(i)**2) * exp( X(1) * params%x_values(i) + X(2) ) &
                  )
             h(2) = &
                  h(2) + f(i)* ( &
                  - params%x_values(i) * exp( X(1) * params%x_values(i) + X(2) ) &
                  )
             h(4) = &
                  h(4) + f(i)* ( &
                  -  exp( X(1) * params%x_values(i) + X(2) ) &
                  )
          end do
          h(3) = h(2)
       end select

       status = 0
       
!!$! let's use Powell's function for now....
!!$       f(1) = X(1) + 10.0 * X(2)
!!$       f(2) = sqrt(5.0) * (X(3) - X(4))
!!$       f(3) = ( X(2) - 2.0 * X(3) )**2
!!$       f(4) = sqrt(10.0) * ( X(1) - X(4) )**2
       
! end of subroutine eval_F
       
     END SUBROUTINE eval_H


     subroutine generate_data_example(x_data,y_data,num_observations)
       
       USE ISO_FORTRAN_ENV
       
       real( real64 ), intent(out) :: x_data(:), y_data(:)
       integer ( int32 ), intent(in) :: num_observations

       ! First, let's get the data
       ! Generated with the code
       !  randn('seed', 23497);
       !   m = 0.3;
       !   c = 0.1;
       !   x_data = [0:0.075:5];
       !   y = exp(m * x_data + c);
       !   noise = randn(size(x_data)) * 0.2;
       !   y_data = y + noise;
       ! (c.f. https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/curve_fitting.cc)
       x_data = (/ 0.0, &
   0.075000000000000, &
   0.150000000000000, &
   0.225000000000000, &
   0.300000000000000, &
   0.375000000000000, &
   0.450000000000000, &
   0.525000000000000, &
   0.600000000000000, &
   0.675000000000000, &
   0.750000000000000, &
   0.825000000000000, &
   0.900000000000000, &
   0.975000000000000, &
   1.050000000000000, &
   1.125000000000000, &
   1.200000000000000, &
   1.275000000000000, &
   1.350000000000000, &
   1.425000000000000, &
   1.500000000000000, &
   1.575000000000000, &
   1.650000000000000, &
   1.725000000000000, &
   1.800000000000000, &
   1.875000000000000, &
   1.950000000000000, &
   2.025000000000000, &
   2.100000000000000, &
   2.175000000000000, &
   2.250000000000000, &
   2.325000000000000, &
   2.400000000000000, &
   2.475000000000000, &
   2.550000000000000, &
   2.625000000000000, &
   2.700000000000000, &
   2.775000000000000, &
   2.850000000000000, &
   2.925000000000000, &
   3.000000000000000, &
   3.075000000000000, &
   3.150000000000000, &
   3.225000000000001, &
   3.300000000000000, &
   3.375000000000000, &
   3.450000000000000, &
   3.525000000000000, &
   3.600000000000001, &
   3.675000000000000, &
   3.750000000000000, &
   3.825000000000000, &
   3.900000000000000, &
   3.975000000000000, &
   4.050000000000001, &
   4.125000000000000, &
   4.200000000000000, &
   4.275000000000000, &
   4.350000000000001, &
   4.425000000000000, &
   4.500000000000000, &
   4.575000000000000, &
   4.650000000000000, &
   4.725000000000001, &
   4.800000000000000, &
   4.875000000000000, &
   4.950000000000000 /)

       y_data = (/ 0.907946872110432, &
   1.199579396036134, &
   1.060092431384317, &
   1.298370500472354, &
   0.952768858414788, &
   1.209665290655204, &
   1.256912538155493, &
   1.163922146095987, &
   1.004877938808100, &
   1.205944250961060, &
   0.952693297695969, &
   1.449662692280761, &
   1.402015259144406, &
   1.378094012325746, &
   1.560882147577552, &
   1.437185539058121, &
   1.559853079888265, &
   1.877814947316832, &
   1.818781749024682, &
   1.375546045112591, &
   1.233967904388409, &
   1.887793124397751, &
   1.610237096463521, &
   1.787032484792262, &
   1.850015127982676, &
   2.120553361509177, &
   1.942913663511919, &
   2.106517132599766, &
   2.271787117356578, &
   1.727554346001754, &
   2.002909500898113, &
   1.975837413903495, &
   2.337446525801909, &
   1.960190841677278, &
   2.447097025572309, &
   2.161663720225506, &
   2.748798529374621, &
   2.507814238594416, &
   2.423769408403069, &
   2.578119353028746, &
   2.460310096221557, &
   2.638362783992324, &
   2.765540456237868, &
   2.837165966564409, &
   3.179711963042789, &
   3.245315453091675, &
   3.289631922410174, &
   3.360995198615834, &
   3.470489725998371, &
   3.169513520153466, &
   3.363740517933189, &
   3.665288099084969, &
   3.620334359722351, &
   4.018911445550667, &
   3.512715166706162, &
   3.874661411575566, &
   4.197746303653517, &
   3.703511523106007, &
   4.076351488309604, &
   4.056340365649961, &
   4.297751562451419, &
   4.373076571153739, &
   4.577093065941748, &
   4.856619059058190, &
   4.927350280596274, &
   4.703122139742729, &
   4.870205182453842 /)
       
     end subroutine generate_data_example

end module example_module
