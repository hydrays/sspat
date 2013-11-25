program main
  use motor
  use random
  use setting
  implicit none
  real dwell_time, dwell_mean, trans_time
  integer i, k
  integer npoint 
  real rmin, rstep
  real tATP, tATP_mean, tATP2

  call read_xdata()

  open(88, file="out/sample.txt")
  open(888, file="out/mean.txt")
  open(99, file="out/trans.txt")

  npoint = 100
  rmin = 0.1
  rstep = 0.1
  do k = 1, npoint
     r = rmin + rstep*(k-1)
     dwell_mean = 0.0
     tATP_mean = 0.0

     do i = 1, nsample
        call unfolding(dwell_time, tATP)
        dwell_mean = dwell_mean + dwell_time
        tATP_mean = tATP_mean + tATP
     end do
     
     call translocation(trans_time, tATP2)

     dwell_mean = dwell_mean / nsample
     tATP_mean = tATP_mean / nsample
     write(888, '(3(f15.2))'), r, dwell_mean, tATP_mean
     write(*, '(3(f15.2))'), r, dwell_mean, tATP_mean

     ! Translocation
     print *, 'translocation...'
     write(99, '(4(f15.2))'), r, trans_time, tATP2, tATP2/trans_time
     write(*, '(4(f15.2))'), r, trans_time, tATP2, tATP2/trans_time
  end do

  close(88)
  close(888)
  close(99)

end program main
