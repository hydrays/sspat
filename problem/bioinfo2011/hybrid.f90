program moran
  use paras
  use random
  use nrtype
  implicit none

  integer x(0:m)
  integer(I4B) i, j, k, is_nag
  real(kind=8) s(0:m)
  real(kind=8) mu(0:m)
  real(kind=8) a(0:m, 0:m)
  integer(I4B) r
  real(kind=8) tau, t, u, sum_s, dte, pt, suma_c, dt
  real(kind=8), parameter :: EPS = 1e-20, epsilon = 0.1
  logical cmask(0:m, 0:m)

  character(len=50) :: filename
  integer(I4B) s_index
  character(len=50) :: forstr1

  write(forstr1, '(A, I2.2, A)'), '(I6, F18.4, ', m+1, 'I12)'
  write(filename, '(A, I1, A, I1, A, I2.2, A, I1)'), &
       'data', int(log10(real(Nsample))), '_hybrid_N', &
       nint(log10(real(N))), 'm', m, &
       'mu', -nint(log10(mu0))
  filename=adjustl(filename)
  filename=trim(filename)
  print *, filename
  open (unit = 11, file=filename, action="write")

  call ran_seed(sequence=1234)

  call makemu(mu)
  call makes(s, delta_s)
  do k = 1, NSample
     t = 0.0
     pt = 0
     x=0
     x(0)=N
     cmask = .true.
     main_loop: do while (x(m).eq.0)
        sum_s = 0
        wait_cevent : do while (.true.)
           call getrate(x, s, mu, a)
           call partition(x, cmask)
           suma_c = sum(a, cmask)
           call expdev(u)
           dte = u / suma_c
           sum_s = sum_s + dte
           if (sum_s > delta_t_coarse) then
              dt = delta_t_coarse
              do i = 0, m
                 do j=0, m
                    if (.not.cmask(i, j)) then
                       r = poidev( a(i, j)*dt )
                       if (i.eq.j) then
                          x(i) = x(i) - r
                          x(i+1) = x(i+1) + r
                          print *, 'mutation as coarse?'
                          read(*, *)
                       else
                          x(i) = x(i) - r
                          x(j) = x(j) + r
                       end if
                    end if
                 end do
              end do
              t = t + dt
              exit wait_cevent
           else 
              call ran1(u)
              u = suma_c*u
              ssa_select: do i = 0, m
                 do j = 0, m
                    if ( cmask(i, j) ) then
                       u = u - a(i, j)
                       if (u < 0) then
                          if (i.eq.j) then
                             x(i) = x(i) - 1
                             x(i+1) = x(i+1) + 1
                          else
                             x(i) = x(i) - 1
                             x(j) = x(j) + 1
                          end if
                          if (x(m).ne.0) then
                             exit main_loop
                          end if
                          exit ssa_select
                       end if
                    end if
                 end do
              end do ssa_select
           end if
        end do wait_cevent
!        if ( k.eq.109 .and. t> 6000) then
!           print *, t, x
!       end if
       ! print *, t, dt, x 
       ! write(11, '(f12.6)'), dt
     end do main_loop
     if (x(m).eq.2) then
        print *, cmask
        print *, t
        read(*,*)
     end if
     write (*, forstr1), k, t, x
     write (11, forstr1), k, t, x
     !read(*, *)
  end do
  close(unit=11)
end program moran
