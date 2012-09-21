module kmcsample3d
  implicit none

contains
  subroutine kmc_sample_omp3d32()
    use omp_lib
    use par_zig_mod
    use setting
    implicit none

    real t, tau, tp, u, private_t
    integer output_index, i, j, active_index, l
    integer k1, k2, shift_i, shift_j 
    integer ibl, jbl
    integer ilow, iup, jlow, jup
    integer nthread, scanner
    real(4) r
    integer :: grainsize = 32
    integer, allocatable :: seed(:)
    real Lbox_npar, Lbox_2npar
    integer ibox, jbox

    Lbox_npar = Lbox/npar
    Lbox_2npar = Lbox/(2*npar)

    print *, Lbox_npar, Lbox_2npar
    !read(*,*)

    allocate(seed(2**npar))
    do i = 1,2**npar
       call random_number(r)
       seed(i) = thisrandseed*r
    enddo
    call par_zigset(2**npar, seed, grainsize)

    open (unit = 100, file='./out/logfile', action="write")

    call init_cell_pool()
  
    t = 0.0
    tp = 0.0
    private_t = 0.0
    output_index = 0
    scanner = 1

    call omp_set_num_threads(int(2**npar, 4))

    do while (t < tend)
       if (t .ge. tp) then
          call cell_stat(t)
          call output_to_file(output_index)
          output_index = output_index + 1
          tp = tp + tpinc
       end if
       t = t + timestep

       if (t .ge. tm) then
          do i = 60, 190
             do j = 60, 190
                do l = 1, npack(i,j)
                   !cmat(i,j,l)%gene3 = 0.0008
                   cmat(i,j,l)%type = 0
                end do
                npack(i,j) = 0
                TDC(i,j) = 0
                call Update_Rate(i, j)
                call Update_Rate(i-1, j)
                call Update_Rate(i+1, j)
                call Update_Rate(i, j-1)
                call Update_Rate(i, j+1)
             end do
             tm = huge(1.0)
          end do
       end if

       do ibl = 1, 0, -1
          do jbl = 1, 0, -1
             !$OMP PARALLEL default(private) &
             !$OMP shared(a, NT, NP, t, npack, TDC, &
             !$OMP cmat, ibl, jbl, scanner, npar, &
             !$OMP Lbox_npar, Lbox_2npar, timestep)

             nthread = OMP_GET_THREAD_NUM()
             ibox = nthread/npar
             jbox = mod(nthread, npar)

             !print *, nthread, ibox, jbox
             !read(*,*)

             ilow = ibox*Lbox_npar + ibl*Lbox_2npar + 1
             iup = ilow + Lbox_2npar - 1

             jlow = jbox*Lbox_npar + jbl*Lbox_2npar + 1
             jup = jlow + Lbox_2npar - 1

             !write(*, '(8I6)'), nthread, ibox, jbox, ilow, iup, jlow, jup
             private_t = t
             !print *, nthread, private_t, t
             !read(*,*)
             do while ( private_t < t +  timestep)
                call Next_Reaction(k1, k2, tau, ilow, iup, jlow, jup)
                !print *, 'next', k1, k2
                if (k1.eq.0 .and. k2.eq.0) then
                   private_t = private_t + 0.00001
                else
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
                end if
             end do
             !$OMP END parallel
             !$OMP barrier
          end do
          !read(*,*)
          !scanner = scanner + 1
       end do
    end do
    close(unit=100)
  end subroutine kmc_sample_omp3d32

  subroutine kmc_sample_omp3d64()
    use omp_lib
    use par_zig_mod
    use setting
    implicit none

    real t, tau, tp, u, private_t
    integer output_index, i, j, active_index, l
    integer k1, k2, shift_i, shift_j 
    integer ibl, jbl
    integer ilow, iup, jlow, jup
    integer nthread, scanner
    real(4) r
    integer :: grainsize = 64
    integer, allocatable :: seed(:)
    real Lbox_npar, Lbox_2npar
    integer ibox, jbox

    Lbox_npar = Lbox/npar
    Lbox_2npar = Lbox/(2*npar)

    !print *, Lbox_npar, Lbox_2npar
    !read(*,*)

    allocate(seed(2**npar))
    do i = 1,2**npar
       call random_number(r)
       seed(i) = thisrandseed*r
    enddo
    call par_zigset(2**npar, seed, grainsize)

    open (unit = 100, file='./out/logfile', action="write")
  
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    private_t = 0.0
    output_index = 0
    scanner = 1

    call omp_set_num_threads(int(2**npar, 4))

    do while (t < tend)
       if (t .ge. tp) then
          call cell_stat(t)
          call output_to_file(output_index)
          output_index = output_index + 1
          tp = tp + 1.0
       end if
       t = t + timestep

       if (t .ge. tm) then
          do i = 60, 190
             do j = 60, 190
                do l = 1, npack(i,j)
                   !cmat(i,j,l)%gene3 = 0.0008
                   cmat(i,j,l)%type = 0
                end do
                npack(i,j) = 0
                TDC(i,j) = 0
                call Update_Rate(i, j)
                call Update_Rate(i-1, j)
                call Update_Rate(i+1, j)
                call Update_Rate(i, j-1)
                call Update_Rate(i, j+1)
             end do
             tm = huge(1.0)
          end do
       end if

       do ibl = 1, 0, -1
          do jbl = 1, 0, -1
             !$OMP PARALLEL default(private) &
             !$OMP shared(a, NT, NP, t, npack, TDC, &
             !$OMP cmat, ibl, jbl, scanner, npar, &
             !$OMP Lbox_npar, Lbox_2npar, timestep)

             nthread = OMP_GET_THREAD_NUM()
             ibox = nthread/npar
             jbox = mod(nthread, npar)

             !print *, nthread, ibox, jbox
             !read(*,*)

             ilow = ibox*Lbox_npar + ibl*Lbox_2npar + 1
             iup = ilow + Lbox_2npar - 1

             jlow = jbox*Lbox_npar + jbl*Lbox_2npar + 1
             jup = jlow + Lbox_2npar - 1

             !write(*, '(8I6)'), nthread, ibox, jbox, ilow, iup, jlow, jup
             private_t = t
             !print *, nthread, private_t, t
             !read(*,*)
             do while ( private_t < t +  timestep)
                call Next_Reaction(k1, k2, tau, ilow, iup, jlow, jup)
                !print *, 'next', k1, k2
                if (k1.eq.0 .and. k2.eq.0) then
                   private_t = private_t + 0.00001
                else
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
                end if
             end do
             !$OMP END parallel
             !$OMP barrier
          end do
          !read(*,*)
          !scanner = scanner + 1
       end do
    end do
    close(unit=100)
  end subroutine kmc_sample_omp3d64

end module kmcsample3d
