module kmcsample
  use setting
  implicit none

contains
  subroutine kmc_sample_serial()
    use random
    use setting
    implicit none

    real t, tau, tp, u, tn
    integer output_index, i, j, active_index
    integer k, shift_i, kill_number
    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    tn = 0.0
    output_index = 0

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          !call cell_stat(t)
          output_index = output_index + 1
          tp = tp + tpinc
       end if

       ! if (t .ge. tm) then
       !    kill_number = 200 
       !    if ( cmat(kill_number, 1)%type .eq. 1 ) then
       !       SC(kill_number) = SC(kill_number) - 1
       !    else if ( cmat(kill_number, 1)%type .eq. 2 ) then
       !       TAC(kill_number) = TAC(kill_number) - 1
       !    else if ( cmat(kill_number, 1)%type .eq. 3 ) then
       !       TDC(kill_number) = TDC(kill_number) - 1
       !    else if ( cmat(kill_number, 1)%type .eq. 5 ) then
       !    else
       !       print *, "do not suppose to find a mutantion cell now...error"
       !       read(*, *)
       !    end if
       !    cmat(kill_number, 1)%type = 4
       !    MC(kill_number) = MC(kill_number) + 1
       !    tm = tm + 10000.0
       ! end if

       ! call Next_Reaction(k, tau)
       
       do i = 1, L
          do j = 1, H
             call cell_event(i, j)
          end do
       end do
       t = t + timestep
       ! if ( (k .le. 2).or.(k .ge. L-1) ) then
       !    call Perodic_BC(k)
       ! end if

       ! do i = 1, L
       !    NT(i) = NT(i) + a(i)*tau
       ! end do

       ! do i = -2, 2
       !    shift_i = k + i
       !    if ( shift_i .le. 0 ) then
       !       shift_i = shift_i + L
       !    else if ( shift_i > L ) then
       !       shift_i = shift_i - L
       !    end if
       !    call Update_Rate(shift_i)
       ! end do

       !call expdev(u)
       ! if ( u.eq.0 ) then
       !    print *, 'exp value = 0'
       !    read(*,*)
       ! end if
       ! NP(k) = NP(k) + u
       ! t = t + tau

    end do
    close(unit=100)

  end subroutine kmc_sample_serial
end module kmcsample
