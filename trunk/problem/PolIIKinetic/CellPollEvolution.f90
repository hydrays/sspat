! This code is abondoned because 
!  1. the cell type is changed. 
!  2. when dividing, the initial status of the daughter cells is changed. 
! in this case, the unbiased enviorment may give biased selection on cells.
! in the modifed case, no biased selection is observed.
! There are some other changes, refer to CellPoolMutation.f90


program ssa
  use setting
  use random
  implicit none

  real u
  real tau, t, tp, td
  integer ClockTime
  integer i, j, k, index, pindex
  real death_toll(NSample)

  call ran_seed(sequence=1234)

  open (unit = 99, file="out/tr.dat", action="write")

  ! Initial value
  do index = 1, NSample
     CellPool(index)%x = xinit
     CellPool(index)%id = mod(index, 3)
  end do
  pindex = 1

  ! Evolution
  do ClockTime = 0, 10000
     call output_to_file2(99, ClockTime)
     if ( mod(ClockTime, 100).eq.0 ) then
        call output_to_file(pindex)
        pindex = pindex + 1
     end if
     ! Update enviornment
     if ( ClockTime > t_now ) then
        if ( env .eq. 0.0 ) then
           env = 10.0
           t_now = t_now + t_on
        else if ( env .eq. 10.0 ) then
           env = 0.0
           t_now = t_now + t_off
        else
           print *, "something wrong"
           read(*,*)
        end if
     end if
     ! Update all cells
     do index = 1, NSample
        call evolve_cell(index, 0.1)      
        if (CellPool(index)%x(5) > v_mature) then
           ! Kill a cell
           print *, "attempting to kill a cell...", index, CellPool(index)%x(5)
           death_toll = 0.0
           death_toll(1) = 1.0 / CellPool(1)%x(5)
           do i = 2, NSample
              death_toll(i) = death_toll(i-1) + 1.0/CellPool(i)%x(5)
           end do
           call ran2(u)
           u = u * death_toll(NSample)
           j = 1
           do while (death_toll(j) .le. u)
              j = j+1
           end do
           print *, "cell j is killed", CellPool(j)%x(5)
           CellPool(index)%x(5) = floor(CellPool(index)%x(5) / 2.0)
           CellPool(j) = CellPool(index)
        end if
     end do
  end do
  close(99)
end program ssa
