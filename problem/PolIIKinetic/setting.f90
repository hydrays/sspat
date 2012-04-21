!!$ Pol II kinetic Simulation

!!$ >>> NSpec = 4

!!$ (1,  2,  3,  4)
!!$  C   R   E  RNA

!!$ >>> NReac = 7

!!$  1: I -> C
!!$  2: C -> I
!!$  3: C -> R
!!$  4: R -> C
!!$  5: R -> E + C
!!$  6: E -> RNA
!!$  7: R -> E

module chem_data  
  use nrtype
  implicit none
  integer(I4B) :: NSample = 1000
  integer(I4B) :: Ntrail = 10000
  integer(I4B), parameter :: NSpec=4
  integer(I4B), parameter :: NReac=7
  real(kind=8) s1, s2, s3
  real(kind=8), parameter :: ep = 0.01
  real NT(NReac)
  real NP(NReac)
  real(kind=8) a(NReac)

  integer(I4B), parameter, dimension(NSpec,NReac) :: nu = reshape( &
       source=(/ &
       (/01, 00, 00, 00/), & !1
       (/-1, 00, 00, 00/), & !2
       (/-1, 01, 00, 00/), & !3
       (/01, -1, 00, 00/), & !4
       (/01, -1, 01, 00/), & !5
       (/00, 00, -1, 01/), & !6
       (/00, -1, 01, 00/) &  !7
       /), shape = (/NSpec, NReac/) &
       )

  real(kind=8), parameter, dimension(NReac) :: c =  &
       (/0.01, 0.0, 1.0, 0.0, 10.0, 0.1, -1.0/)

  real(kind=8), parameter, dimension(NSpec) :: xinit =  &
       (/0.0, 0.0, 0.0, 0.0/)

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
       print *, x
       read(*,*)
    end if

    a(1) = s1*c(1)*is_avai
    a(2) = c(2)*x(1)
    a(3) = s2*c(3)*x(1)
    a(4) = c(4)*x(2)
    a(5) = (1-ep)*s3*c(5)*x(2)
    a(6) = c(6)*x(3)
    a(7) = ep*s3*c(5)*x(2)

  end subroutine getrate

  subroutine Next_Reaction(k, tau)
    implicit none
    integer(I4B), intent(out) :: k
    real(kind=8), intent(out) :: tau
    real(kind=8) tau_temp
    integer(I4B) i

    tau = huge(0.0)
    k = 0
    do i = 1, NReac
       if ( a(i) > 0.0 ) then
          tau_temp = ( NP(i) - NT(i) ) / a(i)
          if ( tau_temp < tau) then
             tau = tau_temp
             k = i
          end if
       else
       end if
    end do
    if ( tau > 0.1 ) then
       tau = 0.1
       k = 0
    end if
    !if ( k .eq. 0 ) then
    !   tau = 0.1
    !end if
    if ( tau < 0 .or. k < 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, k, NP(k) - NT(k), a(k)
       read(*,*)
    end if
  end subroutine Next_Reaction

end module chem_data
