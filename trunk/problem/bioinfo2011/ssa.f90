program ssa
  use paras
  use random
  implicit none

  integer x(0:m)
  integer i, j, counter
  real(kind=8) tau
  real(kind=8) t
  real(kind=8) u
  real(kind=8) s(0:m)
  real(kind=8) mu(0:m)
  real(kind=8) cuma(0:msminusone)
  real(kind=8) a(0:msminusone)
  character(len=50) :: filename
  character(len=50) :: forstr1

  write(forstr1, '(A, I1, A)'), '(I6, F12.4, ', m+1, 'I12)'
  write(filename, '(A, I1, A, I1, A, I1, A, I1, A, I1, A, I1, A, I1)'), &
       'data', int(log10(real(Nsample))), '_ssa_N', int(log10(real(N))), &
       'm', m, 'mu', -nint(log10(mu0))

  filename=adjustl(filename)
  filename=trim(filename)
  print *, filename
  open (unit = 11, file=filename, action="write")	

  call ran_seed(sequence=1234)
  call makemu(mu)
  call makes(s)
  call makenu(nu)
  do i = 1, NSample
     x=0
     x(0)=N
     t=0.0
     counter = 0
     do while (x(m) .eq. 0)
        counter = counter + 1
        call getrate(x, s, mu, a)
        cuma(0) = a(0)
        do j=1, msminusone
           cuma(j) = cuma(j-1) + a(j)
        end do
        call expdev(tau)
        t = t + tau/cuma(msminusone)
        call ran1(u)
        u = u*cuma(msminusone)
        j = 0
        do while ( cuma(j) < u )
           j = j+1
        end do
        x = x + nu(:,j)
     end do
     write (*, forstr1), i, t, x
     write (11, forstr1), i, t, x
  end do
  close(unit=11)
end program ssa
