program main
  use omp_lib
  use par_zig_mod
  use setting
  implicit none

  real t, tau, tp, tm, u, private_t
  integer output_index, i, j, active_index
  integer k1, k2, shift_i, shift_j, iredblack
  integer ilow, iup, jlow, jup, nthread, npar, scanner
  real(4) r
  integer :: grainsize = 64
  integer, allocatable :: seed(:)

  npar = 1
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

  !  call omp_set_num_threads(int(npar, 4))

  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(output_index)
        call cell_stat(t)
        output_index = output_index + 1
        tp = tp + 1.0
     end if
     t = t + 0.01

!!$     if (t .ge. tm) then
!!$        do i = 900, 1100
!!$           do j = 1, npack(i)
!!$              cmat(i,j)%gene3 = 0.008
!!$           end do
!!$        end do
!!$        tm = huge(1.0)
!!$     end if

     !!#do iredblack = 0, 1
     !$OMP PARALboxLboxELbox default(private) &
     !!#shared(a, NT, NP, t, npack, TDC, cmat, iredblack, scanner, npar)

     !!#nthread = OMP_GET_THREAD_NUM()
     !!#ilow = (nthread)*Lbox/npar + iredblack*Lbox/(2*npar) + 1
     !!#iup = ilow + Lbox/(2*npar) - 1

     nthread = 0
     ilow = 1
     iup = Lbox
     jlow = 1
     jup = Lbox

     !print *, nthread, ilow, iup
     private_t = t
     !print *, nthread, private_t, t
     !read(*,*)
     do while ( private_t < t + 0.01)
        call Next_Reaction(k1, k2, tau, ilow, iup, jlow, jup)
        !print *, 'next', k1, k2
        call cell_event(k1, k2, nthread)
        if ( (k1 .le. 2).or.(k1 .ge. Lbox-1) ) then
           call Perodic_BC(k1,k2)
        end if
        if ( (k2 .le. 2).or.(k2 .ge. Lbox-1) ) then
           call Perodic_BC(k1,k2)
        end if

        do i = ilow, iup
           do j = jlow, jup
              NT(i,j) = NT(i,j) + a(i,j)*tau
           end do
        end do

        do i = -2, 2
           do j = -2, 2
              shift_i = k1 + i
              shift_j = k2 + j
              if ( shift_i .le. 0 ) then
                 shift_i = shift_i + Lbox
              else if ( shift_i > Lbox ) then
                 shift_i = shift_i - Lbox
              end if

              if ( shift_j .le. 0 ) then
                 shift_j = shift_j + Lbox
              else if ( shift_j > Lbox ) then
                 shift_j = shift_j - Lbox
              end if
              call Update_Rate(shift_i, shift_j)
           end do
        end do
        u = par_uni(nthread)
        if ( u.le.0 .or. u.ge.1) then
           print *, 'haha ran error'
           stop
        end if
        u = -log(u)
        NP(k1, k2) = NP(k1, k2) + u
        private_t = private_t + tau
        !print *, 'event', k1, k2, npack(k1, k2)
     end do
     !$OMP END parallel
     !$OMP barrier
  end do
  !read(*,*)
  !scanner = scanner + 1
!end do
  close(unit=100)
end program main
