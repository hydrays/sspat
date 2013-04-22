program hybrid
  use paras
  use random
  use nrtype
  implicit none

  integer(I4B) x(0:m)
  integer(I4B) i, j, counter
  real(kind=8) s(0:m)
  real(kind=8) mu(0:m)
  real(kind=8) cuma(0:msminusone)
  real(kind=8) cuma_c(0:msminusone)
  real(kind=8) a(0:msminusone)
  real(kind=8) a_c(0:msminusone)
  real(kind=8) r(0:msminusone)

  real(kind=8) tau
  real(kind=8) t
  real(kind=8) u

  integer(I4B) critical_index(0:msminusone)

  character(len=50) :: filename
  integer(I4B) s_index
  character(len=50) :: forstr1

  write(forstr1, '(A, I2.2, A)'), '(I6, F18.4, ', m+1, 'I12)'
  call ran_seed(sequence=12341)
  write(filename, '(A, I1, A, I1, A, I2.2, A, I1)'), &
       'data', int(log10(real(Nsample))), '_hybrid_N', &
       nint(log10(real(N))), 'm', m, &
       'mu', -nint(log10(mu0))
  filename=adjustl(filename)
  filename=trim(filename)
  print *, filename
  open (unit = 11, file=filename, action="write")	

  call makemu(mu)
  call makes(s, delta_s)
  !  call makes(s)
  call makenu(nu)
  do i = 1, NSample
     x=0
     x(0)=N
     t=0.0
     counter = 0
     do while (x(m) .eq. 0)
        counter = counter +1
        call getrate(x, s, mu, a)
        call partition(x, critical_index)
        call expdev(tau)
        call ran1(u)        
        if (sum(critical_index).eq.ms) then ! only ssa
           cuma(0) = a(0)
           do j=1, msminusone
              cuma(j) = cuma(j-1) + a(j)
           end do
           tau = tau/cuma(msminusone)
           u = u*cuma(msminusone)
           j = 0
           do while ( cuma(j) < u )
              j = j+1
           end do
           x = x + nu(:,j)
        else ! tauleap 
           a_c = a*critical_index
           cuma_c(0) = a_c(0)
           do j=1, msminusone
              cuma_c(j) = cuma_c(j-1) + a_c(j)
           end do
           if(cuma_c(msminusone) .ge. null_a) then ! dirct tauleap
              tau = tau/cuma_c(msminusone)
              u = u*cuma_c(msminusone)
              j = 0
              do while ((critical_index(j).eq.0).or.(cuma_c(j) < u) )
                 j = j+1
              end do
              x = x + nu(:,j)
           else ! no-reaction tauleap
              tau = tau/null_a
              u = u*null_a
              if (u < cuma_c(msminusone)) then
                 j = 0
                 do while ((critical_index(j).eq.0).or.(cuma_c(j) < u) )
                    j = j+1
                 end do
                 x = x + nu(:,j)
              end if
           end if
           r = 0
           do j=0, msminusone
              if(critical_index(j) .eq. 0) then
                 r(j) = poidev( a(j)*tau )
                 x = x + r(j)*nu(:,j)
              end if
           end do
        end if
        t = t + tau
        !write (*, forstr1), i, tau, x
        !write (11, forstr1), i, tau, x
     end do
     !print *, counter
     write (*, forstr1), i, t, x
     write (11, forstr1), i, t, x
     !read(*,*)
  end do
  close(unit=11)
end program hybrid
