program main
  use omp_lib
  use par_zig_mod
  use setting
  implicit none

  real t, tau, tp, tm, u, private_t
  integer output_index, i, j, active_index
  integer k, shift_i, iredblack
  integer ilow, iup, nthread, npar
  real(4) r
  integer :: grainsize = 32
  integer, allocatable :: seed(:)

  npar = 2
  allocate(seed(npar))
  do i = 1,npar
     call random_number(r)
     seed(i) = 123456789*r
  enddo
  call par_zigset(npar, seed, grainsize)

  open (unit = 100, file='./out/logfile', action="write")
  !call ran_seed(sequence=12345)

  call init_cell_pool()

  t = 0.0
  tp = 0.0
  tm = 1000.0
  private_t = 0.0
  output_index = 0
  scanner = 1

  CALL OMP_SET_NUM_THREADS(npar)

  do while (t < tend)
     if (t > tp) then
        call output_to_file(output_index)
        call cell_stat(t)
        output_index = output_index + 1
        tp = tp + 2.0
     end if

     t = t + 0.1

     do iredblack = 0, 1
     !$OMP PARALLEL default(private) &
     shared(a, NT, NP, t, npack, TDC, cmat, iredblack, npar)
     
     nthread = OMP_GET_THREAD_NUM()
     ilow = (nthread)*L/npar + iredblack*L/(2*npar) + 1
     iup = ilow + L/(2*npar) - 1
     
!     print *, nthread, ilow, iup
     private_t = t
     !print *, nthread, private_t, t
     !read(*,*)
     do while ( private_t < t + 0.1)
        call Next_Reaction(k, tau, ilow, iup)
        call cell_event(k, nthread)
        if ( (k .le. 2).or.(k .ge. L-1) ) then
           call Perodic_BC(k)
        end if
        
        do i = ilow, iup
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
        u = par_uni(nthread)
        if ( u.le.0 .or. u.ge.1) then
           print *, 'haha ran error'
           stop
        end if
        u = -log(u)
        NP(k) = NP(k) + u
        private_t = private_t + tau
        !print *, nthread, private_t, 'event', k
     end do
     !$OMP END PARALLEL
     !$OMP barrier
  end do
!  read(*,*)
  end do
  close(unit=100)
end program main
