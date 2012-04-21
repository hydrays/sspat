program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te
  real(kind=8) a(NReac), cuma(NReac), u
  real(kind=8) delta_t, t, tp, td, mtime, ftime
  real(kind=8) takeover_counter, nottakeover_counter, coexist_counter
  real(kind=8) takeover_flag, nottakeover_flag
  real(kind=8) min_TDCnum, average_min_TDCnum
  integer(I4B) is_nag, tac_mflag, sc_mflag
  integer(I4B) i, j, k, index
  real(kind=8) N_mutation
  real(kind=8) xbar(NSpec), tbar, xbar_counter

  call ran_seed(sequence=12341)
  te = 300.0
  pm = 1.0
  vm = 1.0

!  do pm = 0.5, 1.0001, 0.01
!  do vm = 0.1, 3.01, 0.01
  takeover_counter = 0.0
  nottakeover_counter = 0.0
  coexist_counter = 0.0
  average_min_TDCnum = 0.0
  do index = 1.0, NSample
     x = xinit
     t = 0.0
     tp = 0.0
     td = 200.0
     tbar = 0.0
     xbar_counter = 0.0
     takeover_flag = 0.0
     nottakeover_flag = 0.0
     min_TDCnum = huge(0.0)
     ftime = 0.0
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

        call checkx(x, is_nag)
        if (is_nag .eq. 1) then
           print *, 'nag'
           pause
        end if
        
        if(t > tp) then
           write (*, '(F10.2, 10F8.2)'), t, x, sum(x), p0, v0
           tp =  tp + 1.0
        end if

        if(t > td) then
           if ( x(3) .eq. 0 ) then
              x(3) = 10
           end if
           td =  td + 400000.0
        end if

!!$        if ( t > 202 .and. x(2) < min_TDCnum ) then
!!$           min_TDCnum = x(2)
!!$           ftime = t
!!$        end if

!!$        if ( x(1)+x(2).eq.0 ) then
!!$           takeover_counter = takeover_counter + 1.0
!!$           takeover_flag = 1.0
!!$           exit
!!$        end if
!!$
!!$        if ( t > 201.0 .and. x(3)+x(4).eq.0 ) then
!!$           nottakeover_counter = nottakeover_counter + 1.0
!!$           nottakeover_flag = 1.0
!!$           exit
!!$        end if

     end do
     if ( takeover_flag .eq. 0.0 .and. nottakeover_flag .eq. 0.0 ) then
        coexist_counter = coexist_counter + 1.0
     end if
     average_min_TDCnum = average_min_TDCnum + min_TDCnum
!     write (*, '(F10.2, 10F8.2)'), min_TDCnum, ftime
  end do
  average_min_TDCnum = average_min_TDCnum/real(NSample)
!  write (*, '(10F12.2)'), pm, vm, takeover_counter, nottakeover_counter, &
!       coexist_counter, average_min_TDCnum 
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

