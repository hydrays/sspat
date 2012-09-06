!!$ Stem cell simulation

!!$ >>> NSpec = 4
!!$ (1,  2,  3,  4)
!!$  SC  TC  MC  TMC

!!$ >>> NReac = 12
!!$  1: SC -> 2 SC
!!$  2: SC -> SC + TC
!!$  3: SC -> 2 TC
!!$  4: MC -> 2 MC
!!$  5: MC -> MC + TMC
!!$  6: MC -> 2 TMC
!!$  7: TC -> 0
!!$  8: TMC -> 0
!!$  9: SC -> 0
!!$  10: TC -> 0
!!$  11: MC -> 0
!!$  12: TMC -> 0

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 1
  integer(I4B), parameter :: NSpec=4
  integer(I4B), parameter :: NReac=12
  real(kind=8) ap, p0, v0
  real(kind=8) pm, vm
  real(kind=8) ksc1, ksc2, kmc1, kmc2
  real(kind=8) v0max, v0min, vmmax, vmmin
  real(kind=8) q1, q2, q3
  real(kind=8) qm1, qm2, qm3
  real(kind=8), parameter :: L = 200
  real(kind=8), parameter :: mu = 0.0

  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       50, & !SC
       150, & !TC
       00, & !MC
       00 & !TMC
       /)

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       ! 1   2   3   4   5
       source=(/ &
       (/01, 00, 00, 00/), & !1
       (/00, 01, 00, 00/), & !2
       (/-1, 02, 00, 00/), & !3
       (/00, 00, 01, 00/), & !4
       (/00, 00, 00, 01/), & !5
       (/00, 00, -1, 02/), & !6
       (/00, -1, 00, 00/), & !7
       (/00, 00, 00, -1/), & !8
       (/-1, 00, 00, 00/), & !9
       (/00, -1, 00, 00/), & !10
       (/00, 00, -1, 00/), & !11
       (/00, 00, 00, -1/) & !12
       /), shape = (/NSpec, NReac/) &
       )

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    real(kind=8) TGFbeta 

    !TGFbeta = x(2) + x(4)
    TGFbeta = x(2)

    v0max = 3.0
    v0min = 0.5
    ksc1 = 1.0
    ksc2 = v0max/v0min - 1.0
    p0 = 1.0/(1.01 + ksc1*TGFbeta/L)
    v0 = v0max/(1.0 + ksc2*TGFbeta/L)
    !v0 = 1.0!0.65

    !vmmax = 3.0
    !vmmin = 0.5
    !kmc1 = 1.0
    !kmc2 = v0max/v0min - 1.0
    !pm = 1.0/(1.01 + kmc1*TGFbeta/L)
    !vm = vmmax/(1.0 + kmc2*TGFbeta/L)    
    pm = 0.4
    vm = 1.0

!!$    if ( p0 .le. 0.5 ) then
!!$       q2 = 2.0*p0
!!$    else
!!$       q2 = -2.0*(p0 - 1.0)
!!$    end if
!!$    q1 = p0 - q2*0.5
!!$    q3 = 1.0 - q1 - q2
!!$
    q1 = p0
    q2 = 0.0
    q3 = 1.0 - q1
!!$
!!$    if ( pm .le. 0.5 ) then
!!$       qm2 = 2.0*pm
!!$    else
!!$       qm2 = -2.0*(pm - 1.0)
!!$    end if
!!$    qm1 = pm - qm2*0.5
!!$    qm3 = 1.0 - qm1 - qm2
!!$
    qm1 = pm
    qm2 = 0.0
    qm3 = 1.0 - qm1

    if (sum(x) > L) then
       ap = 0.1*(sum(x) - L)
    else
       ap = 0.0
    end if
    a(1) = q1*v0*x(1)
    a(2) = q2*v0*x(1)
    a(3) = q3*v0*x(1)

    a(4) = qm1*vm*x(3)
    a(5) = qm2*vm*x(3)
    a(6) = qm3*vm*x(3)
    
    a(7) = 0.2*x(2)
    a(8) = 0.2*x(4)

    a(9) = ap*x(1)
    a(10) = ap*x(2)
    a(11) = ap*x(3)
    a(12) = ap*x(4)

  end subroutine getrate

end module chem_data
