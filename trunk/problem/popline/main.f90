program main
  use nrtype
  use random
  use setting
  implicit none

  real t, tp, tm, u1
  integer output_index, i
  type(cell), allocatable :: cmat(:,:)

  call ran_seed(sequence=1234)
  allocate(cmat(0:L+1, H))

  cmat(1:L, 1)%type = 1
  cmat(1:L, 2:3)%type = 2
  cmat(1:L, 3:4)%type = 3
  do i = 1, L
     call ran2(u1)
     cmat(i, 1)%gene1 = 0.5*u1
  end do

  t = 0.0
  tp = 0.0
  tm = 1000.0
  output_index = 0
  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(cmat, output_index)
        call cell_stat(cmat, t)
        output_index = output_index + 1
        tp = tp + 10.0
     end if

!!$     if (t .ge. tm) then
!!$        cmat(150:250, :)%type = 0
!!$        tm = tm + 6000.0
!!$     end if

     cmat(0, :) = cmat(L, :)
     cmat(L+1, :) = cmat(1, :)
     call cell_dd(cmat)
     cmat(0, :) = cmat(L, :)
     cmat(L+1, :) = cmat(1, :)
     call cell_mm(cmat)
!!$     cmat(0, :) = cmat(L, :)
!!$     cmat(L+1, :) = cmat(1, :)
!!$     call cell_mutation(cmat, t)
     t = t + delta_t
  end do

  deallocate(cmat)
end program main

subroutine cell_mm(cmat)
  use random
  use setting
  implicit none
  type(cell), intent(inout) :: cmat(0:L+1,H)

  integer i, j, k
  integer pack_num(0:L+1)
  real u1, u2
  real tran_p1, tran_p2
  type(cell) temp


  pack_num = 0
  do i = 1, L
     do j = 1, H
        if ( cmat(i,j)%type .ne. 0 ) then
           pack_num(i) = pack_num(i) + 1
        end if
     end do
     if ( pack_num(i) > 1 ) then
        do k = pack_num(i), 2, -1
           if ( cmat(i,k)%type .eq. 1 ) then
              temp = cmat(i,k-1)
              cmat(i,k-1) = cmat(i,k)
              cmat(i,k) = temp
           end if
        end do
     end if
  end do
  cmat(L+1, :) = cmat(1, :)
  cmat(0, :) = cmat(L, :)

  do i = 1, L
     pack_num = 0
     do j = 1, H
        if ( cmat(i-1,j)%type .ne. 0 ) then
           pack_num(i-1) = pack_num(i-1) + 1
        end if
        if ( cmat(i,j)%type .ne. 0 ) then
           pack_num(i) = pack_num(i) + 1
        end if
        if ( cmat(i+1,j)%type .ne. 0 ) then
           pack_num(i+1) = pack_num(i+1) + 1
        end if
     end do

     if ( pack_num(i) > 1) then
        tran_p1 = D*delta_t * exp(real(pack_num(i)-pack_num(i+1)))
        tran_p2 = D*delta_t * exp(real(pack_num(i)-pack_num(i-1)))
        call ran2(u1)
        if ( u1 < tran_p1 +tran_p2 ) then
           call ran2(u2)
           u2 = u2 * (tran_p1 + tran_p2)
           if ( u2 < tran_p1 ) then
              cmat(i+1,pack_num(i+1)+1)=cmat(i,pack_num(i))
              cmat(i,pack_num(i))%type = 0
           else
              cmat(i-1,pack_num(i-1)+1)=cmat(i,pack_num(i))
              cmat(i,pack_num(i))%type = 0
           end if
        end if
     end if
     
     ! perodic boundary condition
     if ( i.eq. 1 ) then
        cmat(L, :) = cmat(0, :)
        cmat(L+1, :) = cmat(1, :)
     end if
     if ( i.eq. 2 ) then
        cmat(L+1, :) = cmat(1, :)
     end if
     if ( i.eq. L ) then
        cmat(1, :) = cmat(L+1, :)
        cmat(0, :) = cmat(L, :)
     end if
  end do
end subroutine cell_mm

subroutine cell_dd(cmat)
  use random
  use setting
  implicit none

  type(cell), intent(inout) :: cmat(0:L+1,H)
  integer i, j
  real u
  
  call getTGFbeta(cmat)
  
  !$omp parallel private(u)
  do i = 1, L
     do j = 1, H
        call ran2(u)
        if ( u < delta_t*v ) then
           call cell_event(cmat, i, j)
        end if
     end do
  end do
  !$omp end parallel
end subroutine cell_dd

subroutine cell_stat(cmat, t)
  use setting
  implicit none

  type(cell), intent(in) :: cmat(0:L+1,H)
  real, intent(in) :: t

  integer i, j, k
  integer num_sc, num_tac, num_tdc, num_mc

  num_sc = 0
  num_tac = 0
  num_tdc = 0
  num_mc = 0
  do i = 1, L
     do j = 1, H
        if (cmat(i,j)%type.eq.1) then
           num_sc = num_sc + 1
        elseif (cmat(i,j)%type.eq.2) then
           num_tac = num_tac + 1
        elseif (cmat(i,j)%type.eq.3) then
           num_tdc = num_tdc + 1
        elseif (cmat(i,j)%type.eq.4) then
           num_mc = num_mc + 1
        end if
     end do
  end do
  
  write(*, '(5(F10.2))'), t, real(num_sc)/L, &
       real(num_tac)/L, real(num_tdc)/L, real(num_mc)/L
end subroutine cell_stat

!!$subroutine cell_mutation(cmat, t)
!!$  use random
!!$  use setting
!!$  implicit none
!!$  type(cell), intent(inout) :: cmat(0:L+1,H)
!!$  real, intent(in) :: t
!!$
!!$  integer i, j
!!$  real u1, u2
!!$
!!$  call ran2(u1)
!!$  if ( u1 < dt ) then
!!$     print *, 'mutation occured at time ', t
!!$     do i = 1, L
!!$        do j = 1, H
!!$           if (cmat(i,j)%type.eq.1) then
!!$              
!!$           end if
!!$        end do
!!$     end do
!!$  end if
!!$     
!!$end subroutine cell_mutation
