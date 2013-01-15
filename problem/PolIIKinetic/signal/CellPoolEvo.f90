program ssa
  use setting
  use random
  implicit none

  real u
  real tau, t, tp
  integer ClockTime
  integer i, j, k, index, pindex

  call ran_seed(sequence=1234)

  open (unit = 99, file="out/trrb1.dat", action="write")

  ! Initial value
  do index = 1, NSample
     CellPool(index)%x = xinit
     CellPool(index)%id = 1
  end do
  pindex = 0
  tp = 0.0
  exlevel = 0.0
  call update_signal(0.0)

  ! Evolution
  do index = 0, 1200
     t = 0.1*index
     call update_signal(t)
!     if ( 0.1*index .ge. tp ) then
     call output_to_file2(99, t)
     pindex = pindex + 1
     !tp = tp + exwindow
     exlevel = 0.0
     print *, index
     !end if
     ! Update all cells
     !print *, index
     do i = 1, NSample
        call evolve_cell(i, 0.1)      
     end do
  end do
  close(99)
end program ssa
