program main
  use random
  use setting
  implicit none
  real dwell_time, dwell_mean
  real tau
  real nATP
  real jump_prob
  integer jp_flag, is_nag
  real a(2), cuma(2)
  real u, u2
  integer j, i, k
  integer npoint 
  real rmin, rstep
  real tATP, tATP_mean

  call read_xdata()

  if (tr_flag .eq. 1) then
     open(8, file="out/tr.txt")
  end if

  open(88, file="out/sample.txt")
  open(888, file="out/mean.txt")

  npoint = 100
  rmin = 0.02
  rstep = 0.1
  do k = 1, npoint
     r = rmin + rstep*(k-1)
     dwell_mean = 0.0
     tATP_mean = 0.0

     do i = 1, nsample

        jp_flag = 0
        dwell_time = 0.0
        nATP = 0.0
        tATP = 0.0

        do while(jp_flag .eq. 0)
           if (tr_flag .eq. 1) then
              write(8, '(2(f10.2))'), dwell_time, nATP 
              write(*, '(2(f10.2))'), dwell_time, nATP 
           end if
           call getrate(nATP, a)
           cuma = a
           do j=2, 2
              cuma(j) = cuma(j-1) + a(j)
           end do
           call expdev(tau)
           tau = tau/cuma(2)
           dwell_time = dwell_time + tau

           call ran2(u)
           u = cuma(2)*u
           j = 1
           do while (cuma(j) .le. u)
              j = j+1
           end do
           if ( j .eq. 1 ) then
              nATP = nATP + 1.0
           else if ( j .eq. 2 ) then
              ! tempering
              tATP = tATP + nATP
              call ran2(u2)
              if ( u2 < jp_base*exp(nATP*jp_incre) ) then
                 jp_flag = 1.0
                 dwell_mean = dwell_mean + dwell_time
                 tATP_mean = tATP_mean + tATP
                 !write(88, '(1(f10.2))'), dwell_time
                 !write(*, '(1(f10.2))'), dwell_time
              end if
              nATP = 0.0
           else
              print *, 'not implemented'
              read(*,*)
           end if
        end do

     end do
     dwell_mean = dwell_mean / nsample
     tATP_mean = tATP_mean / nsample
     write(888, '(3(f15.2))'), r, dwell_mean, tATP_mean
     write(*, '(3(f15.2))'), r, dwell_mean, tATP_mean
  end do

  if (tr_flag .eq. 1) then
     close(8)
  end if
  close(88)
  close(888)

end program main
