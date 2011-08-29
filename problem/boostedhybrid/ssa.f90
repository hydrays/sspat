program ssa
  use chem_data
  use random
  use nrtype
  implicit none

  real(kind=8) x(NSpec)
  integer(I4B) i, j, k, is_nag
  real(kind=8) cuma(NReac)
  real(kind=8) a(NReac)
  real(kind=8) dt, t, tend, tp
  real(kind=8) u

  tend = end_time
  call ran_seed(sequence=1234)
  do k = 1, NSample
     x = xinit
     t = 0.0
     tp = 1.0
     do while( .true. )
        call getrate(x, a)
        cuma = a
        do j=2, NReac
           cuma(j) = cuma(j-1) + a(j)
        end do
        call expdev(dt)
        t = t + dt/cuma(NReac)
        if (t > tend) then
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
           read(*,*)
        end if
        !if(t > tp) then
        !   write (*, *), t, x
        !   tp =  tp + 1.0
          !read(*,*)
        !end if
     end do
     write(*, '(I6, 10(f8.1))'), k, x
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
