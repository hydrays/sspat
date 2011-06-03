program ssa
  use paras
  use random
  implicit none

  integer i, index
  real(kind=8) tau
  real(kind=8) u1, u2, u3
  real(kind=8) tprint

  type(cell) :: cell_pool(20000) 
  type(cell) :: cell_pool_old(20000) 

  open (unit = 11, file='output', action="write")

  call ran_seed(sequence=1234)
  call pool_init(cell_pool, num)

  t = 0.0
  tprint = 0.0
  index = 0
  do while ( t < tend )
     cell_pool_old = cell_pool
     num_old = num

     next_die_time = cell_pool(1)%cellclock
     next_die_index = 1
     num_cell_can_born = 0
     num_of_tdc = 0
     do i = 1, num
        if ( cell_pool(i)%cellclock < next_die_time ) then
           next_die_time = cell_pool(i)%cellclock
           next_die_index = i
        end if
        if (cell_pool(i)%celltype < 3) then
           num_cell_can_born = num_cell_can_born + 1
        elseif ( cell_pool(i)%celltype .eq. 3 ) then
           num_of_tdc = num_of_tdc + 1
        end if
     end do

     if (num > L1) then
        lambda = lambda0*exp(-0.01*(num - L1))
     else
        lambda = lambda0
     end if
     call expdev(next_born_time)
     next_born_time = next_born_time/(lambda*num_cell_can_born)
     tau = min(next_born_time, next_die_time)
     t =  t + tau
     do i = 1, num
        cell_pool(i)%cellclock = cell_pool(i)%cellclock - tau
     end do
     if ( next_die_time < next_born_time ) then
        cell_pool(next_die_index) = cell_pool(num)
        cell_pool(num)%celltype = 0
        cell_pool(num)%cellclock = 0.
        num = num - 1
        !print *, t, next_die_time
!        read(*,*)
     else
        tau = next_born_time
        call ran2(u1)
        u1 = u1*num_cell_can_born
        next_born_index = 0
        do while ( u1 .ge. 0 )
           next_born_index = next_born_index + 1
           if ( cell_pool(next_born_index)%celltype < 3 ) then
              u1 = u1 - 1
           end if
        end do
        p0 = 1.0/(1.01 + k*real(num_of_tdc)/real(L1))
        p1 = 0.2!1.0/(1.01 + 2.0*k*real(num_of_tdc)/real(L1))
        if ( cell_pool(next_born_index)%celltype .eq. 1 ) then
           call ran2(u3)
           if ( u3 < 0.8 ) then
              p = -1.0
           else
              p = p0
           end if
        elseif ( cell_pool(next_born_index)%celltype .eq. 2 ) then
           !if ( num < L2 ) then
           call ran2(u3)
           if ( u3 < 0.8 ) then
              p = -1.0
           else
              p = p1
           end if
        else
           print *, 'wrong!'
           read(*,*)
        end if
        num = num + 1
        !print *, 'division', p
        if ( p .eq. -1 ) then
           !print *, 'asymm'
           cell_pool(next_born_index)%cellclock = cell_life_length
           cell_pool(num)%celltype = cell_pool(next_born_index)%celltype + 1
           cell_pool(num)%cellclock = cell_life_length      
        else
           call ran2(u2)
           if ( u2 > p ) then
              cell_pool(next_born_index)%celltype = &
                   cell_pool(next_born_index)%celltype + 1
           end if
           cell_pool(next_born_index)%cellclock = cell_life_length
           cell_pool(num) = cell_pool(next_born_index)
           !print *, 'symm division', t, num, next_born_index
           !read(*, *)
        end if
     end if
     if (t > tprint) then
        !call output_to_file(cell_pool_old, num_old, index)
        index = index + 1
        call cell_count(cell_pool, num)
        tprint = tprint + 1.0
     end if
  end do
end program ssa




