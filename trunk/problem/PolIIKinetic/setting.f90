!!$ Pol II kinetic Simulation

!!$ >>> NSpec = 3

!!$ (2,  3,  4)
!!$  C   R   E

!!$ >>> NReac = 5

!!$  1: 0 -> C
!!$  2: C -> R
!!$  3: R -> E
!!$  4: E -> C
!!$  5: E -> 0

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 1
  integer(I4B), parameter :: NSpec=3
  integer(I4B), parameter :: NReac=5
  real(kind=8) s1, s2, s3

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       source=(/ &
       (/01, 00, 00/), & !1
       (/-1, 01, 00/), & !2
       (/00, -1, 01/), & !3
       (/01, 00, -1/), & !4
       (/00, 00, -1/) &  !5
       /), shape = (/NSpec, NReac/) &
       )

  real(kind=8), parameter, dimension(NReac) :: c =  &
       (/0.02, 0.2, 2.0, 0.1, 0.01/)

  real(kind=8), parameter, dimension(NSpec) :: xinit =  &
       (/0.0, 0.0, 0.0/)

contains
  subroutine getrate(x, a)
    implicit none
    real(kind=8), intent(in) :: x(NSpec)
    real(kind=8), intent(out) :: a(NReac)
    real(kind=8) is_avai
    if ( (x(1) + x(2)) .eq. 1.0) then
       is_avai = 0.0
    else if ( (x(1) + x(2)) .eq. 0.0) then
       is_avai = 1.0
    else
       is_avai = -1.0
       print *, 'wrong occupy'
       read(*,*)
    end if

    a(1) = s1*c(1)*is_avai
    a(2) = s2*c(2)*x(1)
    a(3) = s3*c(3)*x(2)
    a(4) = c(4)*x(3)*is_avai
    a(5) = c(5)*x(3)

  end subroutine getrate

end module chem_data
