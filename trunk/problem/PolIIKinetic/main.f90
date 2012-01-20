program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te, u
  real(kind=8) a(NReac), cuma(NReac)
  real(kind=8) delta_t, t, tp, td
  integer(I4B) is_nag
  integer(I4B) i, j, k, index

  call ran_seed(sequence=12341)
  !te = huge(1.0)
  te = 10000.0
  !td = 1000.0
  do index = 1.0, NSample
     x = xinit
     t = 0.0
     tp = 0.0
     s1 = 1.0
     s2 = 1.0
     s3 = 0.001
     do while(.true.)
!!$        if ( t > td ) then 
!!$           !s1 = tiny(1.0)
!!$           !s1 = 1.0
!!$           !s1 = abs(sin(0.005*t))
!!$           !s2 = 1.0
!!$           if ( mod(t, 1000.0) < 1 ) then
!!$              s3 = 1.0
!!$           else
!!$              s3 = 0.0
!!$           end if
!!$        end if
        call getrate(x, a)
        !write(*, '(5(F12.4))'), t, x
        !print *, a
        !read(*,*)
        cuma = a
        do j=2, NReac
           cuma(j) = cuma(j-1) + a(j)
        end do
        if ( cuma(NReac) .eq. 0.0 ) then
           !print *, 'waiting', t
           t = t + 0.1
           cycle
        end if
        call expdev(delta_t)
        delta_t = delta_t/cuma(NReac)
        t = t + delta_t

        if ( t > te ) then
           write(*, '(5(F12.4))'), te, x
           !write(*,*), ' '
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
!!$        if (x(3) .eq. 1) then
!!$           write(*, '(4(F12.4))'), t, x
!!$           exit
!!$        end if
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

