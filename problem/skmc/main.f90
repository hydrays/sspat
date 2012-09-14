program main
  use random
  use setting
  implicit none

  real t, tau, tp, tm, u
  integer output_index, i, j, active_index
  integer k, shift_i
  open (unit = 100, file='./out/logfile', action="write")
  call ran_seed(sequence=12345)

  call init_cell_pool()

  t = 0.0
  tp = 0.0
  tm = 1000.0
  output_index = 0

  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(output_index)
        call cell_stat(t)
        output_index = output_index + 1
        tp = tp + 5.0
     end if

     call Next_Reaction(k, tau)

     call cell_event(k)
     if ( (k .le. 2).or.(k .ge. L-1) ) then
        call Perodic_BC(k)
     end if

     do i = 1, L
        NT(i) = NT(i) + a(i)*tau
     end do

     do i = -2, 2
        shift_i = k + i
        if ( shift_i .le. 0 ) then
           shift_i = shift_i + L
        else if ( shift_i > L ) then
           shift_i = shift_i - L
        end if
        call Update_Rate(shift_i)
     end do

     call expdev(u)
     if ( u.eq.0 ) then
        print *, 'exp value = 0'
        read(*,*)
     end if
     NP(k) = NP(k) + u
     t = t + tau

  end do
  close(unit=100)
end program main
