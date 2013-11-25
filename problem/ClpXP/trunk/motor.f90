module motor
contains

  subroutine unfolding(dwell_time, tATP)
    use random
    use setting
    implicit none
    real, intent(out) :: dwell_time, tATP
    real tau
    real nATP
    real jump_prob
    integer jp_flag, is_nag
    real a(2), cuma(2)
    real u, u2
    integer j

    if (tr_flag .eq. 1) then
       open(8, file="out/tr.txt")
    end if

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
             !write(88, '(1(f10.2))'), dwell_time
             !write(*, '(1(f10.2))'), dwell_time
          end if
          nATP = 0.0
       else
          print *, 'not implemented'
          read(*,*)
       end if
    end do

    if (tr_flag .eq. 1) then
       close(8)
    end if
  end subroutine unfolding

  subroutine translocation(trans_time, tATP)
    use random
    use setting
    implicit none
    real, intent(out) :: trans_time, tATP
    real tau
    real nATP
    integer is_nag
    real a(2), cuma(2)
    real u, u2
    integer j

    if (tr_flag .eq. 1) then
       open(8, file="out/tr.txt")
    end if

    trans_time = 0.0
    nATP = 0.0
    tATP = 0.0

    do while(tATP < 100000)
       if (tr_flag .eq. 1) then
          write(8, '(2(f10.2))'), trans_time, nATP 
          write(*, '(2(f10.2))'), trans_time, nATP 
       end if
       call getrate(nATP, a)
       cuma = a
       do j=2, 2
          cuma(j) = cuma(j-1) + a(j)
       end do
       call expdev(tau)
       tau = tau/cuma(2)
       trans_time = trans_time + tau

       call ran2(u)
       u = cuma(2)*u
       j = 1
       do while (cuma(j) .le. u)
          j = j+1
       end do
       if ( j .eq. 1 ) then
          nATP = nATP + 1.0
       else if ( j .eq. 2 ) then
          tATP = tATP + nATP
          nATP = 0.0
       else
          print *, 'not implemented'
          read(*,*)
       end if
    end do

    if (tr_flag .eq. 1) then
       close(8)
    end if
  end subroutine translocation

end module motor
