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
  integer(I4B) :: NSample = 1
  integer(I4B), parameter :: NSpec=5
  integer(I4B), parameter :: NReac=17
  real(kind=8) ap, p0, p1, v0, q1, q2, q3, vmut
  real(kind=8) k1, k2, v0max, v0min
  real(kind=8) qq1, qq2, qq3
  real(kind=8), parameter :: L = 200
  real(kind=8), parameter :: mu = 0.0

  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       20, & !SC
       30, & !TA
       150, & !TC
       00, & !MC
       00 & !TMC
       /)
!!$
!!$  integer(I4B), parameter :: Xinit(NSpec)=(/ &
!!$       30, & !SC
!!$       50, & !TA
!!$       120, & !TC
!!$       00, & !MC
!!$       00 & !TMC
!!$       /)

!!$  integer(I4B), parameter :: Xinit(NSpec)=(/ &
!!$       100, & !SC
!!$       300, & !TA
!!$       1600, & !TC
!!$       00, & !MC
!!$       00 & !TMC
!!$       /)

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       ! 1   2   3   4   5
       source=(/ &
       (/01, 00, 00, 00, 00/), & !1
       (/00, 01, 00, 00, 00/), & !2
       (/-1, 02, 00, 00, 00/), & !3
       (/00, 01, 00, 00, 00/), & !4
       (/00, 00, 01, 00, 00/), & !5
       (/00, -1, 02, 00, 00/), & !6
!!       (/00, 00, 00, 00, 01/), & !7
!!       (/00, 00, 00, 01, 00/), & !8
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
    v0max = 3.0
    v0min = 0.5
    k1 = 1.0
    k2 = v0max/v0min - 1.0
    p1 = 0.4
    p0 = 1.0/(1.01 + k1*(x(3)+x(5))/L)
    !v0 = v0max/(1.0 + k2*(x(3)+x(5))/L)
    v0 = 0.65

    if ( p0 .le. 0.5 ) then
       q2 = 2.0*p0
    else
       q2 = -2.0*(p0 - 1.0)
    end if
    q1 = p0 - q2*0.5
    q3 = 1.0 - q1 - q2

!!$    q1 = p0
!!$    q2 = 0.0
!!$    q3 = 1.0 - q1

    qq1 = 0.0
    qq2 = 2.0*p1
    qq3 = 1.0 - qq1 - qq2
!!$
!!$    qq1 = p1
!!$    qq2 = 0.0
!!$    qq3 = 1.0 - qq1 - qq2

    if (sum(x) > L) then
       ap = 0.1*(sum(x) - L)
    else
       ap = 0.0
    end if
    a(1) = q1*v0*x(1)
    a(2) = q2*v0*x(1)
    a(3) = q3*v0*x(1)

    a(4) = qq1*x(2)
    a(5) = qq2*x(2)
    a(6) = qq3*x(2)

    a(7) = vmut*x(4)*pm
    a(8) = vmut*x(4)*(1.0-pm)

    !!a(7) = v0*x(5)
    !!a(8) = x(4)

    a(9) = 0.2*x(3)
    a(10) = 0.2*x(5)
    !!a(10) = 0.0*x(5)

    a(11) = ap*x(1)
    a(12) = ap*x(2)
    a(13) = ap*x(3)

    a(14) = ap*x(4)
    a(15) = ap*x(5)

    a(16) = mu*v0*x(1)
    a(17) = mu*x(2)

  end subroutine getrate

end module chem_data
