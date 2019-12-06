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

     ! ! empty
     ! suma = 0.0
     ! do i = 1, Lbox
     !    do j = 1, Lbox
     !       if ( cmat(i,j)%type .eq. 0 ) then
     !          suma = suma + a(i,j)
     !       end if
     !    end do
     ! end do
     ! print *, suma
     ! fire_flag = 0
     ! call random_number(u)
     ! if ( u < suma*dt ) then
     !    u = u/dt
     !    do i = 1, Lbox
     !       do j = 1, Lbox
     !          if ( (fire_flag .eq. 0) .and. (cmat(i,j)%type .eq. 0) ) then
     !             u = u - a(i,j)
     !             if ( u < 0 ) then
     !                adjMcellNum = 0              
     !                do isub = i-R2, i+R2
     !                   if (isub > 0 .and. isub <= Lbox) then
     !                      do jsub = j-R2, j+R2
     !                         if (jsub > 0 .and. jsub <= Lbox) then
     !                            if ( cmat(isub, jsub)%type .eq. 1 ) then
     !                               adjMcellNum = adjMcellNum + 1
     !                            end if
     !                         end if
     !                      end do
     !                   end if
     !                end do
     !                if ( adjMcellNum > 0 ) then
     !                   cmat(i,j)%type = 3
     !                else
     !                   cmat(i,j)%type = 2                    
     !                end if
     !                fire_flag = 1
     !                write(*, *), i, j, adjMcellNum, cmat(i,j)%type
     !             end if
     !          end if
     !       end do
     !    end do
     ! end if

     do i = 1, Lbox
        do j = 1, Lbox
           if ( cmat(i,j)%type == 0 ) then
              call random_number(u)
              if ( u < a(i,j)*dt ) then
                 adjMcellNum = 0              
                 do isub = i-R2, i+R2
                    if (isub > 0 .and. isub <= Lbox) then
                       do jsub = j-R2, j+R2
                          if (jsub > 0 .and. jsub <= Lbox) then
                             if ( cmat(isub, jsub)%type .eq. 1 ) then
                                adjMcellNum = adjMcellNum + 1
                             end if
                          end if
                       end do
                    end if
                 end do
                 if ( adjMcellNum > 0 ) then
                    cmat(i,j)%type = 3
                 else
                    cmat(i,j)%type = 2                    
                 end if
                 !write(*, *), i, j, adjMcellNum, cmat(i,j)%type
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
              end if
           else if ( (cmat(i,j)%type .eq. 2) .or. (cmat(i,j)%type .eq. 3)) then
              ! T cell
              call random_number(u)
              if ( u < beta*dt ) then
                 cmat(i,j)%type = 0
              else
                 call random_number(u)
                 if ( u < moblty*dt ) then
                    phi_temp = -1.0
                    itargt = -1
                    jtargt = -1
                    itest = i-1
                    jtest = j
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                       phi_temp = phi(itest,jtest)
                       itargt = itest
                       jtargt = jtest
                    end if
                    itest = i+1
                    jtest = j
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                       phi_temp = phi(itest,jtest)
                       itargt = itest
                       jtargt = jtest
                    end if
                    itest = i
                    jtest = j-1
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                       phi_temp = phi(itest,jtest)
                       itargt = itest
                       jtargt = jtest
                    end if
                    itest = i
                    jtest = j+1
                    if ( (itest > 0) .and. (itest <= Lbox) .and. (jtest > 0) .and. (jtest <= Lbox) &
                         .and. (cmat(itest, jtest)%type == 0) .and. (phi(itest,jtest) > phi_temp) ) then
                       phi_temp = phi(itest,jtest)
                       itargt = itest
                       jtargt = jtest
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
