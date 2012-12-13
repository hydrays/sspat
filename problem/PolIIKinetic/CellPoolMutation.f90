program ssa
  use setting
  use random
  implicit none

  real u
  real tau, t, tp, td
  integer ClockTime
  integer i, j, k, index, pindex
  real death_toll(NSample)
  real t_now

  call ran_seed(sequence=1234)

  open (unit = 99, file="out/tr.dat", action="write")

  ! Initial value
  pop_ratio = 0.0
  do index = 1, NSample
     CellPool(index)%x = xinit
     CellPool(index)%id = 1+mod(index, 3)
     pop_ratio(CellPool(index)%id) = pop_ratio(CellPool(index)%id) + 1.0
  end do
  pindex = 1
  t_now = t_off
  env = 0.0

  ! Evolution
  do ClockTime = 1, 10000
     call output_to_file2(99, ClockTime)
     if ( mod(ClockTime, 100).eq.0 ) then
        call output_to_file(pindex)
        pindex = pindex + 1
        print *, pindex
     end if
     ! Update enviornment
     call update_env(0.1*real(ClockTime))
     ! Update all cells
     do index = 1, NSample
        call evolve_cell(index, 0.1)      
        if (CellPool(index)%x(5) > v_mature) then
           ! Kill a cell
           !print *, "attempting to kill a cell...", index, CellPool(index)%x(5)
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
           !print *, "cell j is killed", CellPool(j)%x(5)
           pop_ratio(CellPool(j)%id) = pop_ratio(CellPool(j)%id) - 1.0
           pop_ratio(CellPool(index)%id) = pop_ratio(CellPool(index)%id) + 1.0
           CellPool(index)%x(5) = floor(CellPool(index)%x(5) / 2.0)
           CellPool(index)%x(4) = floor(CellPool(index)%x(4) / 2.0)
           CellPool(index)%x(3) = 0.0
           CellPool(j) = CellPool(index)
           if (CellPool(j)%x(1) .eq. 1) then
              CellPool(j)%x(5) = CellPool(j)%x(5) - real(E0)
           else if (CellPool(j)%x(2) .eq. 1) then
              CellPool(j)%x(5) = CellPool(j)%x(5) - 2.0*real(E0)
           end if
           call ran2(u)
           if ( u < 0.25*pmu ) then
              pop_ratio(CellPool(j)%id) = pop_ratio(CellPool(j)%id) - 1.0
              CellPool(j)%id = mod(CellPool(j)%id + 1, 3)+1
              pop_ratio(CellPool(j)%id) = pop_ratio(CellPool(j)%id) + 1.0
           else if ( u < 0.5*pmu ) then
              pop_ratio(CellPool(j)%id) = pop_ratio(CellPool(j)%id) - 1.0
              CellPool(j)%id = mod(CellPool(j)%id + 2, 3)+1
              pop_ratio(CellPool(j)%id) = pop_ratio(CellPool(j)%id) + 1.0
           else if ( u < 0.75*pmu ) then
              pop_ratio(CellPool(index)%id)=pop_ratio(CellPool(index)%id)-1.0
              CellPool(index)%id = mod(CellPool(index)%id + 1, 3)+1
              pop_ratio(CellPool(index)%id)=pop_ratio(CellPool(index)%id)+1.0

           else if ( u .le. pmu ) then
              pop_ratio(CellPool(index)%id)=pop_ratio(CellPool(index)%id)-1.0
              CellPool(index)%id = mod(CellPool(index)%id + 2, 3)+1
              pop_ratio(CellPool(index)%id)=pop_ratio(CellPool(index)%id)+1.0
           !else
              !print *, "not suppose to be here!"
              !read(*,*)
           end if
        end if
     end do
  end do
  close(99)
end program ssa
