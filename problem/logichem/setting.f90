!!$ Stem cell simulation

!!$ >>> NSpec = 5

!!$ (1,  2,  3,  4,  5)
!!$  SC  TA  TC  MC  TMC

!!$ >>> NReac = 8

!!$  1: SC -> 2 SC
!!$  2: SC -> SC + TA
!!$  3: SC -> 2 TA
!!$  4: TA -> 2 TA
!!$  5: TA -> TA + TC
!!$  6: TA -> 2 TC
!!$  7: MC -> 2 MC
!!$  8: MC -> 2 TMC
!!$  9: TC -> 0
!!$  10: TMC -> 0
!!$  11: SC -> 0
!!$  12: TA -> 0
!!$  13: TC -> 0
!!$  14: MC -> 0
!!$  15: TMC -> 0
!!$  16: SC -> MC
!!$  17: TAC -> MC

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 10000
  integer(I4B), parameter :: NSpec=5
  integer(I4B), parameter :: NReac=17
  real(kind=8) ap, p0, p1, v0, symp, symp1
  real(kind=8) kappa, k2
  real(kind=8), parameter :: L = 2000
  real(kind=8), parameter :: mu = 0.000001

  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       100, & !SC
       200, & !TA
       1400, & !TC
       00, & !MC
       00 & !TMC
       /)

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       ! 1   2   3   4   5
       source=(/ &
       (/01, 00, 00, 00, 00/), & !1
       (/00, 01, 00, 00, 00/), & !2
       (/-1, 02, 00, 00, 00/), & !3
       (/00, 01, 00, 00, 00/), & !4
       (/00, 00, 01, 00, 00/), & !5
       (/00, -1, 02, 00, 00/), & !6
       (/00, 00, 00, 01, 00/), & !7
       (/00, 00, 00, -1, 02/), & !8
       (/00, 00, -1, 00, 00/), & !9
       (/00, 00, 00, 00, -1/), & !10
       (/-1, 00, 00, 00, 00/), & !11
       (/00, -1, 00, 00, 00/), & !12
       (/00, 00, -1, 00, 00/), & !13
       (/00, 00, 00, -1, 00/), & !14
       (/00, 00, 00, 00, -1/), & !15
       (/-1, 00, 00, 00, 01/), & !16
       (/00, -1, 00, 01, 00/) & !17
       /), shape = (/NSpec, NReac/) &
       )

contains
  subroutine getrate(x, a, pm)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    real(kind=8), intent(in) :: pm
    kappa = 0.8
    k2 = 4.0
    p1 = 0.4
    p0 = 1.0/(1.01 + kappa*(x(3)+x(5))/L)
    !p1 = 0.5/(1.01 + 0.5*kappa*(x(3)+x(5))/L)
    v0 = 2.5/(1.0 + k2*(x(3)+x(5))/L)
    !p0 = 1.0/(1.01 + kappa*(1450)/L)
    !v0 = 2.5/(1.0 + k2*(1450.0)/L)
    symp = 1.0!/(1.0 + (5.0*kappa*(x(3)+x(5))/L)**2)
    symp1 = 1.0!symp
    if (sum(x) > L) then
       ap = 0.01*(sum(x) - L)
    else
       ap = 0.0
    end if
    a(1) = symp*v0*x(1)*p0
    a(2) = (1.0 - symp)*v0*x(1)
    a(3) = symp*v0*x(1)*(1.0-p0)
    a(4) = symp1*x(2)*p1
    a(5) = (1.0 - symp1)*x(2)
    a(6) = symp1*x(2)*(1.0-p1)

    a(7) = x(4)*pm
    a(8) = x(4)*(1.0-pm)

    a(9) = 0.1*x(3)
    a(10) = 0.1*x(5)

    a(11) = ap*x(1)
    a(12) = ap*x(2)
    a(13) = ap*x(3)

    a(14) = ap*x(4)
    a(15) = ap*x(5)

    a(16) = 0.0*mu*v0*x(1)
    a(17) = 0.0*mu*x(2)

  end subroutine getrate

end module chem_data
