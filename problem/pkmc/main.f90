program main
  use omp_lib
  use nrtype
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
  do i = 1, L
     call Update_Rate(i)
  end do

  do while (t < tend)
!!$     if (t .ge. tm) then
!!$        print *, 'triming'
!!$        do i = 700, 1400
!!$           do j = 1, npack(i)
!!$              if ( cmat(i, j)%type .eq. 3 ) then
!!$                 TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) - D_TGFbeta
!!$              end if
!!$           end do
!!$        end do
!!$        cmat(700:1400, :)%type = 0
!!$        npack(700:1400) = 0
!!$        tm = tm + 600000.0
!!$     end if

     if (t .ge. tp) then
        call output_to_file(output_index)
        call cell_stat(t)
!        print *, t
        output_index = output_index + 1
        tp = tp + 5.0
     end if

     call Next_Reaction(k, tau)

     call cell_event(k)
     if ( (k .le. 2*b+1+1).or.(k .ge. L - 2*b-1) ) then
        call Perodic_BC(k)
     end if

     do i = -2, 2
        shift_i = k + i
        if ( shift_i .le. 0 ) then
           shift_i = shift_i + L
        else if ( shift_i > L ) then
           shift_i = shift_i - L
        end if
        call Update_Rate(shift_i)
     end do

     do i = 1, L
        NT(i) = NT(i) + a(i)*tau
     end do

     call expdev(u)
     NP(k) = NP(k) + u
     t = t + tau

  end do
  close(unit=100)
end program main
