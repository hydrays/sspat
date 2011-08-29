!!$A simple test model for hybrid method.
!!$ system 1 in the CiCP paper
!!$S1 + L1 -> S2 + L2
!!$S2 + L2 -> S1 + L1
!!$S1 -> S3
!!$S3 -> S1
!!$L1 -> 0
!!$L1 -> L2
!!$L2 -> L1

module chem_data  
  use nrtype
  implicit none
  real(kind=8), parameter :: end_time = 100.0
  integer(I4B) :: NSample = 1
  integer(I4B), parameter :: NSpec=5
  integer(I4B), parameter :: NReac=7
  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       1000, & !L1
       1000, & !L2
       1, & !S1
       0, & !S2
       0 & !S3
       /)
  real(kind=8), parameter :: c(NReac)=(/ &
       10.0, & !1
       10.0, & !2
       0.2, & !3
       0.1, & !4
       0.05, & !5
       1.0, & !6
       1.0 & !7
       /)
  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       ! 1   2   3   4   5
       source=(/ &
       (/-1, 01, -1, 01, 00/), & !1
       (/01, -1, 01, -1, 00/), & !2
       (/00, 00, -1, 00, 01/), & !3
       (/00, 00, 01, 00, -1/), & !4
       (/-1, 00, 00, 00, 00/), & !5
       (/-1, 01, 00, 00, 00/), & !6
       (/01, -1, 00, 00, 00/) & !7
       /), shape = (/NSpec, NReac/) &
       )

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    a(1) = c(1)*x(1)*x(3)
    a(2) = c(2)*x(2)*x(4)
    a(3) = c(3)*x(3)
    a(4) = c(4)*x(5)
    a(5) = c(5)*x(1)
    a(6) = c(6)*x(1)
    a(7) = c(7)*x(2)
  end subroutine getrate

end module chem_data
