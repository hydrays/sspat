program main
  use setting
  implicit none

  real t, tau, tp, u
  integer output_index
  integer i, j
  integer isub, jsub
  integer itargt, jtargt, itest, jtest
  real phi_temp
  integer adjTcellNum, adjMcellNum
  real suma
  integer fire_flag
  integer, allocatable :: seed(:)
  integer size
  real roadWeight(5)
  integer moveDirection

  call read_xdata()

  call random_seed(size=size)
  allocate(seed(size))
  seed(:) = iseed
  call random_seed(put=seed)

  open (unit = 100, file='./out/logfile', action="write")

  call init_cell_pool()

  t = 0.0
  tp = 0.0
  output_index = 0

  do while (t < tend)
     if (t .ge. tp) then
        !call cell_stat(t)
        call output_to_file(output_index)
        output_index = output_index + 1
        tp = tp + tpinc
     end if

     do i = 1, Lbox
        do j = 1, Lbox
           if ( cmat(i,j)%type == 0 ) then
              ! empty
              call random_number(u)
              if ( u < a(i,j)*dt ) then
                 cmat(i,j)%type = 2                                     
              end if
           else if ( cmat(i,j)%type .eq. 1 ) then
              ! M cell
              adjTcellNum = 0              
              do isub = i-R1, i+R1
                 if (isub > 0 .and. isub <= Lbox) then
                    do jsub = j-R1, j+R1
                       if (jsub > 0 .and. jsub <= Lbox) then
                          if ( (cmat(isub, jsub)%type .eq. 2).or.(cmat(isub, jsub)%type .eq. 3) ) then
                             adjTcellNum = adjTcellNum + 1
                          end if
                       end if
                    end do
                 end if
              end do
              if ( adjTcellNum > nc ) then
                 cmat(i,j)%type = 0
                 ! activate T cell
                 do isub = i-R1, i+R1
                    if (isub > 0 .and. isub <= Lbox) then
                       do jsub = j-R1, j+R1
                          if (jsub > 0 .and. jsub <= Lbox) then
                             if ( cmat(isub, jsub)%type == 2 ) then
                                cmat(isub, jsub)%type = 3
                             end if
                          end if
                       end do
                    end if
                 end do
              end if
           else if ( cmat(i,j)%type == 3 ) then
              ! activated T cell
              call random_number(u)
              if ( u < beta*dt ) then
                 cmat(i,j)%type = 0
              end if
           else if ( cmat(i,j)%type == 2 ) then
              ! else if ( cmat(i,j)%type == 2 .or. cmat(i,j)%type == 3) then
              ! T cell
              call random_number(u)
              if ( u < beta*dt ) then
                 cmat(i,j)%type = 0
              else
                 call random_number(u)
                 if ( u < moblty*dt ) then
                    ! phi_temp = -1.0
                    ! itargt = -1
                    ! jtargt = -1
                    ! itest = i-1
                    ! jtest = j
                    ! if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                    !      .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                    !    phi_temp = phi(itest,jtest)
                    !    itargt = itest
                    !    jtargt = jtest
                    ! end if
                    ! itest = i+1
                    ! jtest = j
                    ! if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                    !      .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                    !    phi_temp = phi(itest,jtest)
                    !    itargt = itest
                    !    jtargt = jtest
                    ! end if
                    ! itest = i
                    ! jtest = j-1
                    ! if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                    !      .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                    !    phi_temp = phi(itest,jtest)
                    !    itargt = itest
                    !    jtargt = jtest
                    ! end if
                    ! itest = i
                    ! jtest = j+1
                    ! if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                    !      .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                    !    phi_temp = phi(itest,jtest)
                    !    itargt = itest
                    !    jtargt = jtest
                    ! end if

                    ! another method to move the cells
                    roadWeight(:) = 0.0
                    itest = i
                    jtest = j
                    roadWeight(1) = exp(k0*phi(itest,jtest))
                    itest = i-1
                    jtest = j
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) ) then
                       roadWeight(2) = exp(k0*phi(itest,jtest))
                    end if
                    itest = i+1
                    jtest = j
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) ) then
                       roadWeight(3) = exp(k0*phi(itest,jtest))
                    end if
                    itest = i
                    jtest = j-1
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) ) then
                       roadWeight(4) = exp(k0*phi(itest,jtest))
                    end if
                    itest = i
                    jtest = j+1
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) ) then
                       roadWeight(5) = exp(k0*phi(itest,jtest))
                    end if
                    if ( sum(roadWeight) < 1e-8 ) then
                       print *, 'something wrong here'
                    end if
                    call random_number(u)
                    u = u*sum(roadWeight)
                    do moveDirection = 1, 5
                       u = u - roadWeight(moveDirection)
                       if ( u < 0 ) then
                          exit
                       end if
                    end do
                    !print *, t, moveDirection, roadWeight
                    if ( moveDirection == 1 ) then
                       ! does not move
                       itargt = -1
                       jtargt = -1
                    end if
                    if ( moveDirection == 2 ) then
                       itargt = i - 1
                       jtargt = j
                    end if
                    if ( moveDirection == 3 ) then
                       itargt = i + 1
                       jtargt = j
                    end if
                    if ( moveDirection == 4 ) then
                       itargt = i
                       jtargt = j - 1
                    end if
                    if ( moveDirection == 5 ) then
                       itargt = i
                       jtargt = j + 1
                    end if
                    if ( (itargt > 0) .and. (jtargt > 0) ) then
                       ! finish the move
                       cmat(itargt, jtargt)%type = cmat(i,j)%type
                       cmat(i,j)%type = 0
                    end if
                 end if
              end if
           else
              ! not suppose to happen
              !write(*, *), 'ERROR: unknown cell type'
              !stop
           end if
        end do
     end do
     
     call update_phi()
     call update_rate()
     t = t + dt
  end do
  
  close(unit=100)
  
end program main
