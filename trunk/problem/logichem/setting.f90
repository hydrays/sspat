!!$ Stem cell simulation

!!$ >>> NSpec = 5

!!$ (1,  2,  3,  4,  5)
!!$  SC  TA  TC  MC  TMC

!!$ >>> NReac = 18

!!$  1: SC -> 2 SC
!!$  2: SC -> SC + TA
!!$  3: SC -> 2 TA
!!$  4: TA -> 2 TA
!!$  5: TA -> TA + TC
!!$  6: TA -> 2 TC
!!$  7: MC -> 2 MC
!!$  8: MC -> MC + TMC
!!$  9: MC -> 2 TMC
!!$  10: TC -> 0
!!$  11: TMC -> 0
!!$  12: SC -> 0
!!$  13: TA -> 0
!!$  14: TC -> 0
!!$  15: MC -> 0
!!$  16: TMC -> 0
!!$  17: SC -> MC
!!$  18: TAC -> MC

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 1000
  integer(I4B), parameter :: NSpec=5
  integer(I4B), parameter :: NReac=18
  real(kind=8) ap, p0, p1, v0
  real(kind=8) pm, vm
  real(kind=8) k1, k2, v0max, v0min, km1, km2, vmmax, vmmin
  real(kind=8) q1, q2, q3
  real(kind=8) qq1, qq2, qq3
  real(kind=8) qm1, qm2, qm3
  real(kind=8), parameter :: L = 200
  real(kind=8), parameter :: mu = 0.0

!!$  integer(I4B), parameter :: Xinit(NSpec)=(/ &
!!$       10, & !SC
!!$       0, & !TA
!!$       0, & !TC
!!$       00, & !MC
!!$       00 & !TMC
!!$       /)

  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       30, & !SC
       50, & !TA
       120, & !TC
       00, & !MC
       00 & !TMC
       /)

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
       (/00, 00, 00, 01, 00/), & !7
       (/00, 00, 00, 00, 01/), & !8
       (/00, 00, 00, -1, 02/), & !9
       (/00, 00, -1, 00, 00/), & !10
       (/00, 00, 00, 00, -1/), & !11
       (/-1, 00, 00, 00, 00/), & !12
       (/00, -1, 00, 00, 00/), & !13
       (/00, 00, -1, 00, 00/), & !14
       (/00, 00, 00, -1, 00/), & !15
       (/00, 00, 00, 00, -1/), & !16
       (/-1, 00, 00, 01, 00/), & !17
       (/00, -1, 00, 01, 00/) & !18
       /), shape = (/NSpec, NReac/) &
       )

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    real(kind=8) TGFbeta, g1
    
    TGFbeta = x(3) + x(5)
    !TGFbeta = x(3)

    v0max = 3.0
    v0min = 0.5
    !vmmax = 2.0
    vmmin = 0.5
    k1 = 1.0
    !km1 = 0.6
    k2 = v0max/v0min - 1.0
    km2 = vmmax/vmmin - 1.0
    p1 = 0.4

    p0 = 1.0/(1.01 + k1*TGFbeta/L)
    v0 = v0max/(1.0 + k2*TGFbeta/L)
    !v0 = 0.65

    pm = 1.0/(1.01 + km1*TGFbeta/L)
    vm = vmmax/(1.0 + km2*TGFbeta/L)
    !vm = 0.85

    !v0 = 0.65
    !v0 = 1.0
!!$

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

    if ( pm .le. 0.5 ) then
       qm2 = 2.0*pm
    else
       qm2 = -2.0*(pm - 1.0)
    end if
    qm1 = pm - qm2*0.5
    qm3 = 1.0 - qm1 - qm2

!!$    qm1 = pm
!!$    qm2 = 0.0
!!$    qm3 = 1.0 - qm1

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

    a(7) = qm1*vm*x(4)
    a(8) = qm2*vm*x(4)
    a(9) = qm3*vm*x(4)
    
    a(10) = 0.2*x(3)
    a(11) = 0.2*x(5)

    a(12) = ap*x(1)
    a(13) = ap*x(2)
    a(14) = ap*x(3)

    a(15) = ap*x(4)
    a(16) = ap*x(5)

    a(17) = mu*v0*x(1)
    a(18) = mu*x(2)

  end subroutine getrate

end module chem_data
