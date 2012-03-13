program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te
  real(kind=8) a(NReac), cuma(NReac), u
  real(kind=8) delta_t, t, tp, td, tac_mtime, sc_mtime
  integer(I4B) is_nag, tac_mflag, sc_mflag
  integer(I4B) i, j, k, index
  real(kind=8) pm, N_mutation
  real(kind=8) xbar(NSpec), tbar, xbar_counter

  call ran_seed(sequence=1234)
  te = 200.0
  !te = huge(1.0)
  pm = 1.0
  vmut = 1.0

  do pm = 0.0, 1.0001, 0.02
!  do vmut = 0.1, 3.1, 0.1
!  xbar = 0.0
  do index = 1.0, NSample
     x = xinit
     t = 0.0
     tp = 0.0
     td = 100.0
     tbar = 0.0
     tac_mflag = 0
     sc_mflag = 0
     tac_mtime = 0.0
     sc_mtime = 0.0
     xbar_counter = 0.0
     do while(t < te)
        call getrate(x, a, pm)
        cuma = a
        do j=2, NReac
           cuma(j) = cuma(j-1) + a(j)
        end do
        call expdev(delta_t)
        delta_t = delta_t/cuma(NReac)
        t = t + delta_t
        call ran2(u)
        u = cuma(NReac)*u
        j = 1
        do while (cuma(j) .le. u)
           j = j+1
        end do
        x = x + nu(:, j)
!!$        if (j.eq.16) then
!!$           print *, x
!!$           read(*,*)
!!$        end if
        call checkx(x, is_nag)
        if (is_nag .eq. 1) then
           print *, 'nag'
           pause
        end if
        
!!$        if ( j .eq. 17 .and. tac_mflag.eq.0) then
!!$           tac_mtime = t
!!$           tac_mflag = 1
!!$        end if
!!$
!!$        if ( j .eq. 16 .and. sc_mflag .eq. 0) then
!!$           sc_mtime = t
!!$           sc_mflag = 1
!!$           exit
!!$        end if
!!$        if(t > td) then
!!$           if (x(4).eq.0) x(4) = 1
!!$           td =  td + 200.0
!!$        end if
        if(t > tp) then
           write (*, '(F10.2, 10F8.2)'), t, x, sum(x), ap, p0, v0
           tp =  tp + 1.0
           if (t > 1000) then
           end if
        end if
        !if (x(1).eq.0 .or. x(5).ne.0 .or. x(4).ne.0) then 
        !if (x(1).eq.0 .or. x(5).ne.0) then 
!!$        if (sum(x(1:3)).eq.0) then 
!!$           write (*, '(F10.2, 10F8.2)'), t, x, sum(x), ap, p0, v0           
!!$           exit
!!$        end if
     end do
     xbar = xbar + x
!     write (*, '(F10.2, 9F8.2, 2I8, 2F10.2)'), t, x, sum(x), ap, p0, v0, &
!          tac_mflag, sc_mflag, tac_mtime, sc_mtime
  end do
!  xbar = xbar / Nsample
!  write (*, '(F18.8, 8F10.2)'), pm, vmut, xbar, sum(xbar)
!  end do 
!  end do
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

