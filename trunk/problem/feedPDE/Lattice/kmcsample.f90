module kmcsample
  use setting
  implicit none

contains
  subroutine kmc_sample_serial()
    use random
    use setting
    implicit none

    real t, tau, tp, u
    integer output_index, i, j, active_index
    integer k, shift_i
    real t_update_nutri
    integer kill_number
    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    t_update_nutri = 0.0
    output_index = 0

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          call cell_stat(t)
          output_index = output_index + 1
          tp = tp + tpinc
       end if

       if (t .ge. t_update_nutri .and. t .ge. 50.0 ) then
          call update_nutri(NutriTimestep)
          t_update_nutri = t_update_nutri + NutriTimestep
       end if

       if (t .ge. tm) then
          kill_number = 500 
          if ( cmat(kill_number, 1)%type .eq. 1 ) then
             SC(kill_number) = SC(kill_number) - 1
          else if ( cmat(kill_number, 1)%type .eq. 2 ) then
             TAC(kill_number) = TAC(kill_number) - 1
          else if ( cmat(kill_number, 1)%type .eq. 3 ) then
             TDC(kill_number) = TDC(kill_number) - 1
          else
             print *, "do not suppose to find a mutantion cell now...error"
             read(*, *)
          end if
          cmat(kill_number, 1)%type = 4
          MC(kill_number) = MC(kill_number) + 1
          tm = tm + 10000.0
       end if

       call Next_Reaction(k, tau)

       call cell_event(k)
       if ( (k .le. 2).or.(k .ge. L-1) ) then
          call Perodic_BC(k)
       end if
       ! Check
       if ( MC(k) .ne. 0 .and. k < 300) then
          print *, "error at", k
          print *, TDC(k-3:k+3)
          print *, SC(k-3:k+3)
          print *, TAC(k-3:k+3)
          print *, MC(k-3:k+3)
          print *, npack(k-3:k+3)
          read(*,*)
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

  end subroutine kmc_sample_serial

  subroutine kmc_sample_omp32()
    use omp_lib
    use par_zig_mod
    use setting
    implicit none

    real t, tau, tp, u, private_t
    integer output_index, i, j, active_index
    integer k, shift_i, iredblack
    integer ilow, iup, nthread
    real(4) r
    integer :: grainsize = 32
    integer, allocatable :: seed(:)

    allocate(seed(npar))
    do i = 1,npar
       call random_number(r)
       seed(i) = iseed*r
    enddo
    call par_zigset(npar, seed, grainsize)

    open (unit = 100, file='./out/logfile', action="write")

    call init_cell_pool()

    t = 0.0
    tp = 0.0
    private_t = 0.0
    output_index = 0

    CALL OMP_SET_NUM_THREADS(npar)

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          call cell_stat(t)
          output_index = output_index + 1
          tp = tp + tpinc
       end if

       if (t .ge. tm) then
          i = 500
          j = 1
          cmat(i,j)%type = 4
          cmat(i,j)%gene1 = 0.0
          call Update_Rate(i)
          call Update_Rate(i+1)
          call Update_Rate(i-1)
          tm = huge(1.0)
       end if

       do iredblack = 0, 1
          !$OMP PARALLEL default(private) &
          !$OMP shared(a, NT, NP, t, npack, TDC, &
          !$OMP cmat, iredblack, npar, timestep, L, H)

          nthread = OMP_GET_THREAD_NUM()
          ilow = (nthread)*L/npar + iredblack*L/(2*npar) + 1
          iup = ilow + L/(2*npar) - 1

          !print *, nthread, ilow, iup
          private_t = t
          !print *, nthread, timestep, private_t, t
          !read(*,*)
          do while ( private_t < t + timestep)
             call Next_Reaction_omp(k, tau, ilow, iup)
             call cell_event_omp(k, nthread)
             if ( (k .le. 2).or.(k .ge. L-1) ) then
                call Perodic_BC(k)
             end if
             
             !print *, k, nthread
             
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
       !read(*,*)
       t = t + timestep
    end do
    close(unit=100)

  end subroutine kmc_sample_omp32

  subroutine kmc_sample_omp64()

    use omp_lib
    use par_zig_mod
    use setting
    implicit none

    real t, tau, tp, u, private_t
    integer output_index, i, j, active_index
    integer k, shift_i, iredblack
    integer ilow, iup, nthread
    real(4) r
    integer :: grainsize = 64
    integer, allocatable :: seed(:)

    allocate(seed(npar))
    do i = 1,npar
       call random_number(r)
       seed(i) = iseed*r
    enddo
    call par_zigset(npar, seed, grainsize)

    open (unit = 100, file='./out/logfile', action="write")

    call init_cell_pool()

    t = 0.0
    tp = 0.0
    private_t = 0.0
    output_index = 0

    CALL OMP_SET_NUM_THREADS(npar)

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          call cell_stat(t)
          output_index = output_index + 1
          tp = tp + tpinc
       end if

       if (t .ge. tm) then
          i = 500
          j = 1
          cmat(i,j)%type = 4
          cmat(i,j)%gene1 = 0.0
          call Update_Rate(i)
          call Update_Rate(i+1)
          call Update_Rate(i-1)
          tm = huge(1.0)
       end if

       do iredblack = 0, 1
          !$OMP PARALLEL default(private) &
          !$OMP shared(a, NT, NP, t, npack, TDC, &
          !$OMP cmat, iredblack, npar, timestep, L, H)

          nthread = OMP_GET_THREAD_NUM()
          ilow = (nthread)*L/npar + iredblack*L/(2*npar) + 1
          iup = ilow + L/(2*npar) - 1

          !print *, nthread, ilow, iup
          private_t = t
          !print *, nthread, private_t, t
          !read(*,*)
          do while ( private_t < t + timestep)
             call Next_Reaction_omp(k, tau, ilow, iup)
             call cell_event_omp(k, nthread)
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
       !read(*,*)
       t = t + timestep
    end do
    close(unit=100)

  end subroutine kmc_sample_omp64

end module kmcsample
