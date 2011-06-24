program main
  use nrtype
  use random
  use setting
  implicit none

  real t, tp, tm, u1
  integer output_index, i
  type(cell), allocatable :: cmat(:,:,:)

  call ran_seed(sequence=1234)
  allocate(cmat(0:L+1,0:L+1,H))

  cmat(1:L, 1:L, 1)%type = 1
  cmat(1:L, 1:L, 2:3)%type = 2
  cmat(1:L, 1:L, 3:4)%type = 3
!!$  do i = 1, L
!!$     call ran2(u1)
!!$     cmat(i, 1)%gene1 = 0.5*u1
!!$  end do

  t = 0.0
  tp = 0.0
  tm = 1000.0
  output_index = 0
  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(cmat, output_index)
        call cell_stat(cmat, t)
        output_index = output_index + 1
        tp = tp + 1.0
     end if

!!$     if (t .ge. tm) then
!!$        cmat(150:250, :)%type = 0
!!$        tm = tm + 6000.0
!!$     end if

     cmat(1:L, 0, :) = cmat(1:L, L, :)
     cmat(1:L, L+1, :) = cmat(1:L, 1, :)
     cmat(0, 1:L, :) = cmat(L, 1:L, :) 
     cmat(L+1, 1:L, :) = cmat(1, 1:L, :)
     call cell_dd(cmat)
     cmat(1:L, 0, :) = cmat(1:L, L, :)
     cmat(1:L, L+1, :) = cmat(1:L, 1, :)
     cmat(0, 1:L, :) = cmat(L, 1:L, :) 
     cmat(L+1, 1:L, :) = cmat(1, 1:L, :)
     call cell_mm(cmat)
!!$     cmat(0, :) = cmat(L, :)
!!$     cmat(L+1, :) = cmat(1, :)
!!$     call cell_mutation(cmat, t)
     t = t + delta_t
     print *, t
     !call cell_stat(cmat, t)
  end do

  deallocate(cmat)
end program main

subroutine cell_mm(cmat)
  use random
  use setting
  implicit none
  type(cell), intent(inout) :: cmat(0:L+1,0:L+1,H)

  integer ix, iy, j, k
  integer pack_num(0:L+1, 0:L+1)
  real u1, u2
  real tran_p(4)
  type(cell) temp

  pack_num = 0
  do ix = 1, L
     do iy = 1, L
        do j = 1, H
           if ( cmat(ix,iy,j)%type .ne. 0 ) then
              pack_num(ix,iy) = pack_num(ix,iy) + 1
           end if
        end do
        if ( pack_num(ix,iy) > 1 ) then
           do k = pack_num(ix,iy), 2, -1
              if ( cmat(ix,iy,k)%type .eq. 1 ) then
                 temp = cmat(ix,iy,k-1)
                 cmat(ix,iy,k-1) = cmat(ix,iy,k)
                 cmat(ix,iy,k) = temp
              end if
           end do
        end if
     end do
  end do
  cmat(1:L, 0, :) = cmat(1:L, L, :)
  cmat(1:L, L+1, :) = cmat(1:L, 1, :)
  cmat(0, 1:L, :) = cmat(L, 1:L, :) 
  cmat(L+1, 1:L, :) = cmat(1, 1:L, :)

  do ix = 1, L
     do iy = 1, L
        pack_num = 0
        do j = 1, H
           if ( cmat(ix-1,iy,j)%type .ne. 0 ) then
              pack_num(ix-1,iy) = pack_num(ix-1,iy) + 1
           end if
           if ( cmat(ix+1,iy,j)%type .ne. 0 ) then
              pack_num(ix+1,iy) = pack_num(ix+1,iy) + 1
           end if
           if ( cmat(ix,iy,j)%type .ne. 0 ) then
              pack_num(ix,iy) = pack_num(ix,iy) + 1
           end if
           if ( cmat(ix,iy-1,j)%type .ne. 0 ) then
              pack_num(ix,iy-1) = pack_num(ix,iy-1) + 1
           end if
           if ( cmat(ix,iy+1,j)%type .ne. 0 ) then
              pack_num(ix,iy+1) = pack_num(ix,iy+1) + 1
           end if
     end do

     if ( pack_num(ix,iy) > 1) then
        tran_p(1) = D*delta_t * exp(real(pack_num(ix,iy)-pack_num(ix+1,iy)))
        tran_p(2) = D*delta_t * exp(real(pack_num(ix,iy)-pack_num(ix-1,iy)))
        tran_p(3) = D*delta_t * exp(real(pack_num(ix,iy)-pack_num(ix,iy+1)))
        tran_p(4) = D*delta_t * exp(real(pack_num(ix,iy)-pack_num(ix,iy-1)))
        call ran2(u1)
        if ( u1 < sum(tran_p) ) then
           call ran2(u2)
           u2 = u2 * sum(tran_p)
           u2 = u2 - tran_p(1)
           if ( u2 < 0 ) then ! move right
              cmat(ix+1,iy,pack_num(ix+1,iy)+1)=cmat(ix,iy,pack_num(ix,iy))
              cmat(ix,iy,pack_num(ix,iy))%type = 0
           else 
              u2 = u2 - tran_p(2)
              if ( u2 < 0 ) then ! move left
                 cmat(ix-1,iy,pack_num(ix-1,iy)+1)=cmat(ix,iy,pack_num(ix,iy))
                 cmat(ix,iy,pack_num(ix,iy))%type = 0
              else
                 u2 = u2 - tran_p(3)
                 if ( u2 < 0 ) then ! move up
                    cmat(ix,iy+1,pack_num(ix,iy+1)+1) = &
                         cmat(ix,iy,pack_num(ix,iy))
                    cmat(ix,iy,pack_num(ix,iy))%type = 0
                 else ! move down
                    cmat(ix,iy-1,pack_num(ix,iy-1)+1) = &
                         cmat(ix,iy,pack_num(ix,iy))
                    cmat(ix,iy,pack_num(ix,iy))%type = 0                    
                 end if
              end if
           end if
        end if
     end if
     
     ! perodic boundary condition
     if ( ix .eq. 1 ) then
        cmat(L, iy, :) = cmat(0, iy, :)
        cmat(L+1, iy, :) = cmat(1, iy, :)
     end if
     if ( ix .eq. 2 ) then
        cmat(L+1, iy, :) = cmat(1, iy, :)
     end if
     if ( ix .eq. L ) then
        cmat(1, iy, :) = cmat(L+1, iy, :)
        cmat(0, iy, :) = cmat(L, iy, :)
     end if

     if ( iy .eq. 1 ) then
        cmat(ix, L, :) = cmat(ix, 0, :)
        cmat(ix, L+1, :) = cmat(ix, 1, :)
     end if
     if ( iy .eq. 2 ) then
        cmat(ix, L+1, :) = cmat(ix, 1, :)
     end if
     if ( iy .eq. L ) then
        cmat(ix, 1, :) = cmat(ix, L+1, :)
        cmat(ix, 0, :) = cmat(ix, L, :)
     end if
  end do
end do
end subroutine cell_mm

subroutine cell_dd(cmat)
  use random
  use setting
  implicit none

  type(cell), intent(inout) :: cmat(0:L+1,0:L+1,H)
  integer ix, iy, j
  real u
  
  call getTGFbeta(cmat)
  
  !$omp parallel private(u)
  do ix = 1, L
     do iy = 1, L
        do j = 1, H
           call ran2(u)
           if ( u < delta_t*v ) then
              call cell_event(cmat, ix, iy, j)
           end if
        end do
     end do
  end do
  !$omp end parallel
end subroutine cell_dd

subroutine cell_stat(cmat, t)
  use setting
  implicit none

  type(cell), intent(in) :: cmat(0:L+1,0:L+1,H)
  real, intent(in) :: t

  integer ix,iy, j, k
  integer num_sc, num_tac, num_tdc, num_mc

  num_sc = 0
  num_tac = 0
  num_tdc = 0
  num_mc = 0
  do ix = 1, L
     do iy = 1, L
        do j = 1, H
           if (cmat(ix,iy,j)%type.eq.1) then
              num_sc = num_sc + 1
           elseif (cmat(ix,iy,j)%type.eq.2) then
              num_tac = num_tac + 1
           elseif (cmat(ix,iy,j)%type.eq.3) then
              num_tdc = num_tdc + 1
           elseif (cmat(ix,iy,j)%type.eq.4) then
              num_mc = num_mc + 1
           end if
        end do
     end do
  end do
  write(*, '(5(F10.2))'), t, real(num_sc)/(L*L), &
       real(num_tac)/(L*L), real(num_tdc)/(L*L), real(num_mc)/(L*L)
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
