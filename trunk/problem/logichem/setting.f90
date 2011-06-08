!!$ Stem cell simulation

!!$ >>> NSpec = 3

!!$ (1,  2,  3)
!!$  SC  TA  TC

!!$ >>> NReac = 8

!!$  1: SC -> 2 SC
!!$  2: SC -> 2 TA
!!$  3: TA -> 2 TA
!!$  4: TA -> 2 TC
!!$  5: TC -> 0
!!$  6: SC -> 0
!!$  7: TA -> 0
!!$  8: TC -> 0

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 1
  integer(I4B), parameter :: NSpec=3
  integer(I4B), parameter :: NReac=8
  real(kind=8) ap, p0, p1, v0
  real(kind=8) kappa, k2
  real(kind=8), parameter :: L = 2000

  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       10, & !SC
       00, & !TA
       00 & !TC
       /)

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       ! 1   2   3
       source=(/ &
       (/01, 00, 00/), & !1
       (/-1, 02, 00/), & !2
       (/00, 01, 00/), & !3
       (/00, -1, 02/), & !4
       (/00, 00, -1/), & !5
       (/-1, 00, 00/), & !6
       (/00, -1, 00/), & !7
       (/00, 00, -1/) & !8
       /), shape = (/NSpec, NReac/) &
       )

  real(kind=8), parameter :: end_time = 20.0

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    kappa = 1.0
    k2 = 8.0
    p1 = 0.45
    p0 = 1.0/(1.01+kappa*x(3)/L)
    v0 = 2.0/(1.0 + k2*x(3)/L)
    if (sum(x) > L) then
       ap = 0.001*(sum(x) - L)
    else
       ap = 0.0
    end if
    a(1) = v0*x(1)*p0
    a(2) = v0*x(1)*(1.0-p0)
    a(3) = x(2)*p1
    a(4) = x(2)*(1.0-p1)
    a(5) = 0.2*x(3)
    a(6) = ap*x(1)
    a(7) = ap*x(2)
    a(8) = ap*x(3)
  end subroutine getrate

end module chem_data
