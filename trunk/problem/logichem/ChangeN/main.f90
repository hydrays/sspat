program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te
  real(kind=8) a(NReac), cuma(NReac), u
  real(kind=8) delta_t, t, tp, td, mtime
  integer(I4B) is_nag, tac_mflag, sc_mflag
  integer(I4B) i, j, k, index
  real(kind=8) xbar(NSpec), tbar, xbar_counter
  integer(I4B) FID1, CeilN, flag1
  real(kind=8) dt_wound, dm_wound
  
!  real(kind=8), parameter :: L_space(19) = (/40, 60, 80, 100, &
!       120, 140, 160, 180, 200, 220, 240, 260, 280, &
!       300, 320, 340, 360, 380, 400/)

  real(kind=8), parameter :: L_space(7) = (/40, 60, 100, 200, &
       1000, 2000, 20000/)

!, 3000, 4000, 6000, 10000, &
!       20000, 40000, 80000, 100000, 200000/)


!  real(kind=8), parameter :: L_space(3) = (/800, 1000, 1200/)

  call ran_seed(sequence=12345)

  te = 4000.0
  pm = 1.0
  vm = 1.0

  Fid1 = 10
  dt_wound = 10000.0
  dm_wound = 0.5

  do CeilN = 7,7
     L = L_space(CeilN)
     Xinit = (/2, 3, 15, 0, 0/)
     !Xinit = Xinit*(10**(CeilN-1))
     Xinit = Xinit*L/20.0
     !print *, CeilN, L, Xinit
     call prepare_output_files("mute_pernew_N", L, Fid1, dt_wound, dm_wound)
     do index = 1.0, NSample
        x = xinit
        t = 0.0
        tp = 0.0
        td = 150.0
        flag1 = 0
        tbar = 0.0
        xbar = 0.0
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

              call checkx(x, is_nag)
              if (is_nag .eq. 1) then
                 print *, 'nag'
                 pause
              end if

              if(t > tp) then
                 write (*, '(F10.2, 10F10.2)'), t, x
                 !write (Fid1, '(F10.2, 10F10.2)'), t, x, sum(x), p0, v0
                 tp =  tp + 1.0
              end if

           if(t > td) then
              if ( x(4) .eq. 0 .and. sum(x).ne.0 .and. flag1.eq.0) then
                 x(4) = 1
                 flag1 = 1
                 !x(3) = int(x(3)/4.0)
                 !x(2) = int(x(2)/4.0)
                 !x(1) = int(x(1)/4.0)
              end if
              x(3) = int(dm_wound*x(3))
              td =  td + dt_wound
           end if

           if (sum(x).eq.0) then
              exit
           end if
        end do
        !write (*, '(I10, 10F10.2)'), index, x
        write (Fid1, '(I10, 10F10.2)'), index, x
     end do
     close(Fid1)
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

