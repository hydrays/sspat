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

  call ran_seed(sequence=12341)
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
  integer(I4B), intent(inout) :: x(NSpec)
  integer(I4B),  intent(out) :: is_nag
  is_nag = 0
  if(any(x < 0) ) then
     is_nag = 1
     print *, x
     print *, 'nag!'
  end if
end subroutine checkx

