program ssa
  use paras
  use random
  implicit none

  integer i, index, itype
  real(kind=8) tau
  real(kind=8) u1, u2, u3
  real(kind=8) tprint
  real(kind=8) lambda

  type(cell) :: pool(20000) 
  type(cell) :: pool_old(20000) 

  open (unit = 11, file='output', action="write")

  call ran_seed(sequence=1234)
  call pool_init(pool, num)

  t = 0.0
  tprint = 0.0
  index = 0

  do while ( t < tend )
     pool_old = pool
     num_old = num

     die_time = pool(1)%clock
     die_index = 1
     num_cell_can_born = 0
     num_of_tdc = 0
     do i = 1, num
        if ( pool(i)%clock < die_time ) then
           die_time = pool(i)%clock
           die_index = i
        end if
        if (pool(i)%type < 3) then
           num_cell_can_born = num_cell_can_born + 1
        elseif ( pool(i)%type .eq. 3 ) then
           num_of_tdc = num_of_tdc + 1
        end if
     end do

     call get_lambda(num, lambda)

     call expdev(born_time)
     born_time = born_time/(lambda*num_cell_can_born)

     tau = min(born_time, die_time)
     t =  t + tau
     do i = 1, num
        pool(i)%clock = pool(i)%clock - tau
     end do

     if ( die_time < born_time ) then
        pool(die_index) = pool(num)
        pool(num)%type = 0
        pool(num)%clock = 0.
        num = num - 1
     else
        call ran2(u1)
        u1 = u1*num_cell_can_born
        born_index = 0
        do while ( u1 .ge. 0 )
           born_index = born_index + 1
           if ( pool(born_index)%type < 3 ) then
              u1 = u1 - 1
           end if
        end do

        call getp

        num = num + 1
        
        itype = pool(born_index)%type
        call ran2(u3)
        if ( u3 < 0.0 ) then
           pool(born_index)%clock = lifespan(itype)
           pool(num)%type = pool(born_index)%type + 1
           pool(num)%clock = lifespan(itype+1)      
        else
           call ran2(u2)
           if ( u2 > p(itype) ) then
              itype = itype + 1
              pool(born_index)%type = itype
           end if
           pool(born_index)%clock = lifespan(itype)
           pool(num) = pool(born_index)
        end if
     end if
     if (t > tprint) then
        !call output_to_file(pool_old, num_old, index)
        index = index + 1
        call cell_count(pool, num)
        tprint = tprint + 1.0
     end if
  end do
end program ssa




