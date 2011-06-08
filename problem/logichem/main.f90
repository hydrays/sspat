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
  integer(I4B) i, j, k

  call ran_seed(sequence=1234)

  te = 2000
  do i = 1, Nsample
     x = xinit
     t = 0.0
     tp = 1.0
     td = 0.0
     do while(.true.)
        call getrate(x, a)
        cuma = a
        do j=2, NReac
           cuma(j) = cuma(j-1) + a(j)
        end do
        call expdev(delta_t)
        t = t + delta_t/cuma(NReac)
        if (t > te) then
           exit
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
        if(t > tp) then
           write (*, '(F10.2, 7F10.2, 2E10.2)'), t, x, p0, sum(x), ap, v0
!           print *, 'birth', sum(a(1:4))
!           print *, 'death', sum(a(5:8))
           tp =  tp + 1.0
           !read(*,*)
        end if
        if(t > td) then
           if (x(4).eq.0) x(4) = 1
           td =  td + 200.0
        !        read(*,*)
        end if
     end do
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

