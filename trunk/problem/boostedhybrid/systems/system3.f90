!!$ A simple gene expression model

!!$ >>> NSpec = 7

!!$ (1,  2,  3, 4, 5, 6,  7)
!!$  D   D*  M  R  P  P2  Q)

!!$ >>> NReac = 10

!!$  1: D* -> D + M + R
!!$  2: M -> M + P
!!$  3: M -> 0
!!$  4: P -> 0
!!$  5: D + R -> D*
!!$  6: D* -> D + R
!!$  7: P + P -> P2
!!$  8: P2 -> P + P
!!$  9: D + P2 -> Q
!!$ 10: Q -> D + P2

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 100000
  integer(I4B), parameter :: NSpec=7
  integer(I4B), parameter :: NReac=10
  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       0, & ! D
       1, & ! D*
       0, & ! M
       30, & ! R
       0, & ! P
       0, & ! P2
       0 &  ! Q
       /)
  real(kind=8), parameter :: c(NReac)=(/ &
       0.0078, & !1
       0.043, & !2
       0.0039, & !3
       0.0007, & !4
       0.38, & !5
       3.0, & !6
       0.5, & !7
       5.0, & !8
       0.12, & !9
       9.0 & !10
       /)
  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
!!$ (1,  2,  3, 4, 5, 6,  7)
!!$  D   D*  M  R  P  P2  Q)
       source=(/ &
       ! 1   2   3   4   5   6   7 
       (/01, -1, 01, 01, 00, 00, 00/), & !1
       (/00, 00, 00, 00, 01, 00, 00/), & !2
       (/00, 00, -1, 00, 00, 00, 00/), & !3
       (/00, 00, 00, 00, -1, 00, 00/), & !4
       (/-1, 01, 00, -1, 00, 00, 00/), & !5
       (/01, -1, 00, 01, 00, 00, 00/), & !6
       (/00, 00, 00, 00, -2, 01, 00/), & !7
       (/00, 00, 00, 00, 02, -1, 00/), & !8
       (/-1, 00, 00, 00, 00, -1, 01/), & !9
       (/01, 00, 00, 00, 00, 01, -1/) & !10
       /), shape = (/NSpec, NReac/) &
       )

  real(kind=8), parameter :: end_time = 50000.0

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    a(1) = c(1)*x(2)
    a(2) = c(2)*x(3)
    a(3) = c(3)*x(3)
    a(4) = c(4)*x(5)
    a(5) = c(5)*x(1)*x(4)
    a(6) = c(6)*x(2)
    a(7) = c(7)*x(5)*(x(5)-1)/2.0
    a(8) = c(8)*x(6)
    a(9) = c(9)*x(1)*x(6)
    a(10) = c(10)*x(7)
  end subroutine getrate

end module chem_data
