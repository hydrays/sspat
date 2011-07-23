program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te
  real(kind=8) a(NReac), cuma(NReac), u
  real(kind=8) delta_t, t, tp, td
  integer(I4B) is_nag
  integer(I4B) i, j, k, index
  real(kind=8) pm, N_mutation
  real(kind=8) xbar(NSpec), tbar

  call ran_seed(sequence=1234)
  te = 1000.0
  pm = 1.0
  do index = 1.0, NSample
     N_mutation = 0.0
     x = xinit
     t = 0.0
     tp = 1.0
     td = 50.0
     xbar = 0.0
     tbar = 0.0
     do while(t < te)
        call getrate(x, a, pm)
        cuma = a
        do j=2, NReac
           cuma(j) = cuma(j-1) + a(j)
        end do
        call expdev(delta_t)
        delta_t = delta_t/cuma(NReac)
        t = t + delta_t
        if ( t > 1000 ) then
           xbar = xbar + x*delta_t
           tbar = tbar + delta_t
        end if
        call ran2(u)
        u = cuma(NReac)*u
        j = 1
        do while (cuma(j) .le. u)
           j = j+1
        end do
        x = x + nu(:, j)
        call checkx(x, is_nag)
        if (is_nag .eq. 1) then
           print *, 'nag'
           pause
        end if
        if(t > td) then
           if (x(4).eq.0) x(3) = 50
           td =  td + 100.0
           !read(*,*)
        end if
        if(t > tp) then
           write (*, '(F10.2, 7F10.2, 4E10.2)'), t, x, p0, sum(x), ap, v0, symp
           tp =  tp + 1.0
!           if (t > 400) read(*,*)
        end if
     end do
     xbar = xbar / tbar
     !if (x(4).eq.0) then
     !   N_mutation = N_mutation + 1
     !end if
     !write (*, '(F10.2, 7F10.2)'), t, x, sum(x), N_mutation
     !write (*, '(F18.8, 6F10.2)'), pm, xbar, sum(xbar)
     !read(*,*)
  end do
end program ssa

subroutine checkx(x, is_nag)
  use chem_data
  use nrtype
  implicit none
<<<<<<< .working
  integer(I4B), intent(inout) :: x(NSpec)
  integer(I4B),  intent(out) :: is_nag
  is_nag = 0
  if(any(x < 0) ) then
     is_nag = 1
     print *, x
     print *, 'nag!'
  end if
end subroutine checkx
=======
  integer(I4B) :: NSample = 1000
  integer(I4B), parameter :: NSpec=5
  integer(I4B), parameter :: NReac=17
  real(kind=8) ap, p0, p1, v0, symp, symp1
  real(kind=8) k1, k2, k3
  real(kind=8), parameter :: L = 500
  real(kind=8), parameter :: mu = 0.00001
>>>>>>> .merge-right.r32

<<<<<<< .working
=======
  integer(I4B), parameter :: Xinit(NSpec)=(/ &
       50, & !SC
       100, & !TA
       350, & !TC
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
    k1 = 0.65
    k2 = 4.0
    k3 = 3.0
    p1 = 0.4
    p0 = 1.0/(1.01 + k1*(x(3)+x(5))/L)
    !v0 = 1.0
    v0 = 2.5/(1.0 + k2*(x(3)+x(5))/L)
    !p0 = 1.0/(1.01 + kappa*(1450)/L)
    !v0 = 2.5/(1.0 + k2*(1450.0)/L)
    !p0 = 0.99*exp(-kappa*(x(3)+x(5))/L)
    symp = 1.0/(1.0 + (k3*(x(3)+x(5))/L)**2)
    symp1 = 1.0!symp
    if (sum(x) > L) then
       ap = 0.04*(sum(x) - L)
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

    a(9) = 0.2*x(3)
    a(10) = 0.2*x(5)

    a(11) = ap*x(1)
    a(12) = ap*x(2)
    a(13) = ap*x(3)

    a(14) = ap*x(4)
    a(15) = ap*x(5)

    a(16) = mu*v0*x(1)
    a(17) = mu*x(2)

  end subroutine getrate

end module chem_data
>>>>>>> .merge-right.r32
