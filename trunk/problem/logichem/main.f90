program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te
  real(kind=8) a(NReac), cuma(NReac), u, u1, u2
  real(kind=8) delta_t, t, tp, td, tac_mtime, sc_mtime
  integer(I4B) is_nag, tac_mflag, sc_mflag
  integer(I4B) i, j, k, index
  real(kind=8) N_mutation
  real(kind=8) xbar(NSpec), tbar, xbar_counter

  call ran_seed(sequence=12341)
  te = 1000.0
  !te = huge(1.0)
  !pm = 1.0
  !vmut = 1.0

!  do pm = 0.0, 1.0001, 0.02
!  do vmut = 0.1, 3.1, 0.1
!  xbar = 0.0
  do index = 1.0, NSample
     x = 0
     x(1) = 20
     x(2) = 30
     x(3) = 150

     pm = 0.0
     vmut = 0.0
     
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
        call getrate(x, a)
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
        if (j.eq.41) then
           call ran2(u1)
           call ran2(u2)
           u1 = 0.5 + (1.0-0.5)*u1
           u2 = 0.5 + (2.0-0.5)*u2
!           print *, 'mutation', u1, u2
!           read(*,*)
           if ( x(4) .eq. 0 ) then
              pm(1) = u1
              vmut(1) = u2
              x(4) = 1
           else if ( x(6) .eq. 0 ) then
              pm(2) = u1
              vmut(2) = u2
              x(6) = 1
           else if ( x(8) .eq. 0 ) then
              pm(3) = u1
              vmut(3) = u2
              x(8) = 1
           else if ( x(10) .eq. 0 ) then
              pm(4) = u1
              vmut(4) = u2
              x(10) = 1
           else if ( x(12) .eq. 0 ) then
              pm(5) = u1
              vmut(5) = u2
              x(12) = 1
           else if ( x(14) .eq. 0 ) then
              pm(6) = u1
              vmut(6) = u2
              x(14) = 1
           else
              print *, 'too many mcs...'
              read(*,*)
           end if
        end if
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
!!$!           if (x(4).eq.0) x(4) = 1
!!$           x(3) = 20
!!$           td =  td + 200.0
!!$        end if
!!$        if(t > tp) then
!!$           write (*, '(F10.2, 20F8.2)'), t, x, sum(x), ap, p0, v0
!!$           tp =  tp + 1.0
!!$           if (t > 1000) then
!!$           end if
!!$        end if
        !if (x(1).eq.0 .or. x(5).ne.0 .or. x(4).ne.0) then 
        !if (x(1).eq.0 .or. x(5).ne.0) then 
        if (sum(x(1:3)).eq.0) then 
           write (*, '(F10.2, 20F10.4)'), t, x, sum(x), ap, &
                (pm(1)*x(4)+pm(2)*x(6)+pm(3)*x(8)+pm(4)*x(10)+pm(5)*x(12)+pm(6)*x(14))/(x(4)+x(6)+x(8)+x(10)+x(12)+x(14)), &
                (vmut(1)*x(4)+vmut(2)*x(6)+vmut(3)*x(8)+vmut(4)*x(10)+vmut(5)*x(12)+vmut(6)*x(14))/(x(4)+x(6)+x(8)+x(10)+x(12)+x(14))
           exit
        end if
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

