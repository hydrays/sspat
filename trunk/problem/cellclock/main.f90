program ssa
  use paras
  use random
  implicit none

  integer i, index, itype
  real(kind=8) tau
  real(kind=8) u1, u2, u3
  real(kind=8) tprint, tmutation
  real(kind=8) lambda
  real(kind=8) t

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
     num_of_dc = 0
     do i = 1, num
        if ( pool(i)%clock < die_time ) then
           die_time = pool(i)%clock
           die_index = i
        end if
        if ( (pool(i)%type.eq.3) .or. (pool(i)%type.eq.5) ) then
           num_of_dc = num_of_dc + 1
        end if
     end do

     call getp
     call getv
     call get_lambda(num, lambda)

     born_ruler = 0.0
     do i = 1, num
        itype = pool(i)%type
        if ( (itype.eq.1).or.(itype.eq.2).or.(itype.eq.4) ) then
           born_ruler = born_ruler + v(itype)
        end if
     end do
     call expdev(born_time)
     born_time = born_time/(lambda*born_ruler)
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
        u1 = u1*born_ruler
        born_index = 0
        do while ( u1 .ge. 0 )
           born_index = born_index + 1
           itype = pool(born_index)%type
        if ( (itype.eq.1).or.(itype.eq.2).or.(itype.eq.4) ) then
              u1 = u1 - v(itype)
           end if
        end do

        num = num + 1        
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
        !call output_to_file(pool_old, num_old, index, t)
        index = index + 1
        call cell_count(pool, num, t)
        tprint = tprint + 0.5
     end if
     if (t > tmutation) then
        pool(num)%type = 4
        pool(num)%clock = lifespan(4)
        tmutation = tmutation + 200
     end if
  end do
end program ssa




