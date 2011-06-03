program main
  use nrtype
  use random
  use setting
  implicit none

  real t, tp, tm
  integer output_index
  integer, allocatable :: cell_matrix(:,:)

  call ran_seed(sequence=1234)
  allocate(cell_matrix(0:L+1, D))

  ! 0 --- empty
  ! 1 --- stem cell
  ! 2 --- TAC
  ! 3 --- TDC
  ! 4 --- MC
  ! 5 --- DMC
  cell_matrix = 0
  cell_matrix(1:L, 1) = 1
  cell_matrix(1:L, 2:3) = 2
  cell_matrix(1:L, 3:4) = 3

  t = 0.0
  tp = 0.0
  tm = 500.0
  output_index = 0
  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(cell_matrix, output_index)
        call cell_stat(cell_matrix, t)
        output_index = output_index + 1
        tp = tp + 5.0
!        write(*,'(50(I1))'), cell_matrix(0,:)
!        write(*,'(50(I1))'), cell_matrix(L,:)
!        print *, t, '==========================================='
!        write(*,'(50(I1))'), cell_matrix(1,:)
!        write(*,'(50(I1))'), cell_matrix(2,:)
     end if
     if (t .ge. tm) then
        cell_matrix(150:200, 1:10) = 4
        tm = tm + 6000.0
     end if
     cell_matrix(0, :) = cell_matrix(L, :)
     cell_matrix(L+1, :) = cell_matrix(1, :)
     call cell_dd(cell_matrix)
     cell_matrix(0, :) = cell_matrix(L, :)
     cell_matrix(L+1, :) = cell_matrix(1, :)
     call cell_mm(cell_matrix)
     t = t + delta_t
  end do

  deallocate(cell_matrix)
end program main

subroutine cell_mm(cell_matrix)
  use random
  use setting
  implicit none
  integer, intent(inout) :: cell_matrix(0:L+1,D)

  integer i, j, k
  integer pack_num(0:L+1)
  real u1, u2
  real tran_p1, tran_p2
  integer temp


  pack_num = 0
  do i = 1, L
     do j = 1, D
        if ( cell_matrix(i,j) .ne. 0 ) then
           pack_num(i) = pack_num(i) + 1
        end if
     end do
!!$     if ( pack_num(i) > 1 ) then
!!$        do k = pack_num(i), 2, -1
!!$           if ( cell_matrix(i,k) .eq. 1 .or. cell_matrix(i,k).eq.4 ) then
!!$              temp = cell_matrix(i,k-1)
!!$              cell_matrix(i,k-1) = cell_matrix(i,k)
!!$              cell_matrix(i,k) = temp
!!$           end if
!!$        end do
!!$     end if
  end do
  cell_matrix(L+1, :) = cell_matrix(1, :)
  cell_matrix(0, :) = cell_matrix(L, :)

  do i = 1, L

     pack_num = 0

     do j = 1, D
        if ( cell_matrix(i-1,j) .ne. 0 ) then
           pack_num(i-1) = pack_num(i-1) + 1
        end if
        if ( cell_matrix(i,j) .ne. 0 ) then
           pack_num(i) = pack_num(i) + 1
        end if
        if ( cell_matrix(i+1,j) .ne. 0 ) then
           pack_num(i+1) = pack_num(i+1) + 1
        end if
     end do

     if ( pack_num(i) > 1) then
        if ( cell_matrix(i, pack_num(i)).ne.4 ) then
           tran_p1 = vdiff_high*delta_t * exp(real(pack_num(i)-pack_num(i+1)))
           tran_p2 = vdiff_high*delta_t * exp(real(pack_num(i)-pack_num(i-1)))
           call ran2(u1)
           if ( u1 < tran_p1 +tran_p2 ) then
              call ran2(u2)
              u2 = u2 * (tran_p1 + tran_p2)
              if ( u2 < tran_p1 ) then
                 cell_matrix(i+1,pack_num(i+1)+1)=cell_matrix(i,pack_num(i))
                 cell_matrix(i,pack_num(i)) = 0
              else
                 cell_matrix(i-1,pack_num(i-1)+1)=cell_matrix(i,pack_num(i))
                 cell_matrix(i,pack_num(i)) = 0
              end if
           end if
        else
           tran_p1 = vdiff_low*delta_t * exp(real(pack_num(i)-pack_num(i+1)))
           tran_p2 = vdiff_low*delta_t * exp(real(pack_num(i)-pack_num(i-1)))
           call ran2(u1)
           if ( u1 < tran_p1 +tran_p2 ) then
              call ran2(u2)
              u2 = u2 * (tran_p1 + tran_p2)
              if ( u2 < tran_p1 ) then
                 cell_matrix(i+1,pack_num(i+1)+1)=cell_matrix(i,pack_num(i))
                 cell_matrix(i,pack_num(i)) = 0
              else
                 cell_matrix(i-1,pack_num(i-1)+1)=cell_matrix(i,pack_num(i))
                 cell_matrix(i,pack_num(i)) = 0
              end if
           end if
        end if
     end if
     
     ! perodic boundary condition
     if ( i.eq. 1 ) then
        cell_matrix(L, :) = cell_matrix(0, :)
        cell_matrix(L+1, :) = cell_matrix(1, :)
     end if
     if ( i.eq. 2 ) then
        cell_matrix(L+1, :) = cell_matrix(1, :)
     end if
     if ( i.eq. L ) then
        cell_matrix(1, :) = cell_matrix(L+1, :)
        cell_matrix(0, :) = cell_matrix(L, :)
     end if
  end do
end subroutine cell_mm

subroutine cell_dd(cell_matrix)
  use random
  use setting
  implicit none

  integer, intent(inout) :: cell_matrix(0:L+1,D)

  integer i, j, k
  real temp_num

  real u1, u2, u3

  p0 = 0.
  TGFbeta = 0.
  ! get TGFbeta
  do i = 1, L
     temp_num = 0.
     do j = 1, D
        if ( cell_matrix(i,j) .eq. 3 ) then
           temp_num = temp_num + 1.
        end if
     end do
     if (temp_num > 0) then
        do k = 0, b-1
           if ( k .eq. 0) then
              TGFbeta(i) = TGFbeta(i) + temp_num
           else
              TGFbeta(i+k) = TGFbeta(i+k) + temp_num*(1.0 - real(k)/b)
              TGFbeta(i-k) = TGFbeta(i-k) + temp_num*(1.0 - real(k)/b)
           end if
        end do
     end if
  end do
  do k = 1, b
     TGFbeta(k) = TGFbeta(k) + TGFbeta(L+k)
     TGFbeta(L-k+1) = TGFbeta(L-k+1) + TGFbeta(1-k)
  end do

  do i = 1, L
     p0(i) = 1.0 / (1.0 + 0.008*TGFbeta(i))
  end do

  do i = 1, L
     v0(i) = v0m / (1.0 + 0.1*TGFbeta(i))
  end do

  p1 = 0.4

  do i = 1, L
     do j = 1, D
        if (cell_matrix(i,j) .eq. 1 ) then
           call ran2(u1)
           if ( u1 < delta_t*v0(i) ) then
              ! division
              do k=D, j+2, -1
                 cell_matrix(i, k) = cell_matrix(i, k-1)
              end do
              call ran2(u2)
              if ( u2 < p0(i) ) then
                 ! SC -> 2SC
                 cell_matrix(i, j+1) = 1
              else
                 ! SC -> 2TAC
                 cell_matrix(i, j) = 2
                 cell_matrix(i, j+1) = 2
              end if
           end if
        else if ( cell_matrix(i,j) .eq. 2 ) then
           call ran2(u1)
           if ( u1 < delta_t*v_high ) then
              ! division
              do k=D, j+2, -1
                 cell_matrix(i, k) = cell_matrix(i, k-1)
              end do
              call ran2(u2)
              if ( u2 < p1(i) ) then
                 ! TAC -> 2TAC
                 cell_matrix(i, j+1) = 2
              else
                 ! SC -> 2TAC
                 cell_matrix(i, j) = 3
                 cell_matrix(i, j+1) = 3
              end if
           end if
        else if ( cell_matrix(i,j) .eq. 3 ) then
           call ran2(u1)
           if ( u1 < delta_t*v_high ) then
              ! death
              do k=j, D-1
                 cell_matrix(i, k) = cell_matrix(i, k+1)
              end do
           end if
        else if ( cell_matrix(i,j) .eq. 4 ) then
           call ran2(u1)
           if ( u1 < 0.1*delta_t*v_high ) then
              call ran2(u2)
              if ( u2 < 0.4 ) then
                 ! death
                 do k=j, D-1
                    cell_matrix(i, k) = cell_matrix(i, k+1)
                 end do
              else
                 ! MC division
                 do k=D, j+2, -1
                    cell_matrix(i, k) = cell_matrix(i, k-1)
                 end do
                 call ran2(u3)
                 if ( .true. ) then
                    ! MC -> 2MC
                    cell_matrix(i, j+1) = 4
                 else
                    ! do nothing
                 end if
              end if
           end if
        else
           ! do nothing
        end if
     end do
  end do

end subroutine cell_dd

subroutine cell_stat(cell_matrix, t)
  use setting
  implicit none

  integer, intent(in) :: cell_matrix(0:L+1,D)
  real, intent(in) :: t

  integer i, j, k
  integer num_sc, num_tac, num_tdc, num_mc

  num_sc = 0
  num_tac = 0
  num_tdc = 0
  num_mc = 0
  do i = 1, L
     do j = 1, D
        if (cell_matrix(i,j).eq.1) then
           num_sc = num_sc + 1
        elseif (cell_matrix(i,j).eq.2) then
           num_tac = num_tac + 1
        elseif (cell_matrix(i,j).eq.3) then
           num_tdc = num_tdc + 1
        elseif (cell_matrix(i,j).eq.4) then
           num_mc = num_mc + 1
        end if
     end do
  end do
  
  write(*, '(5(F10.2))'), t, real(num_sc)/L, &
       real(num_tac)/L, real(num_tdc)/L, real(num_mc)/L
end subroutine cell_stat
