module setting
  integer, parameter :: L = 60
  integer, parameter :: H = 60
  real, parameter :: b = 8.0
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 200.0
  real, parameter :: p1 = 0.3
  real, parameter :: v = 1.0
  real, parameter :: D = 1.0
  type cell
     integer type
     real gene1
  end type cell
  type(cell) cmat(0:L+1,H)
  real a(1:L)
  real NT(1:L)
  real NP(1:L)
  integer npack(0:L+1)
  real TGFbeta(-b:L+b+1)
  real D_TGFbeta(1:1+2*b)
  real vr, vl

contains
  subroutine init_cell_pool()
    use random
    implicit none
    integer temp_num
    real u
    integer i, j

    cmat(1:L, 1)%type = 1
    cmat(1:L, 2:3)%type = 2
    cmat(1:L, 3:4)%type = 3
    cmat(0, :) = cmat(L, :)
    cmat(L+1, :) = cmat(1, :)

    npack = 0
    do i = 1, L
       do j = 1, H
          if ( cmat(i,j)%type .ne. 0 ) then
             npack(i) = npack(i) + 1
          end if
       end do
    end do
    npack(0) = npack(L)
    npack(L+1) = npack(1)

    NP = 0.0
    NT = 0.0
    do i = 1, L
       call expdev(u)
       NP(i) = u
    end do

    D_TGFbeta = 0.0
    do j = 1, 1+2*b
       D_TGFbeta(j) = exp(-real(abs(j-b-1))/b)
    end do

    TGFbeta = 0.
    do i = 1, L
       temp_num = 0.
       do j = 1, H
          if ( cmat(i,j)%type .eq. 3 ) then
             TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) + D_TGFbeta
          end if
       end do
    end do
    do j = 1, b
       TGFbeta(j) = TGFbeta(j) + TGFbeta(L+j)
       TGFbeta(L-j+1) = TGFbeta(L-j+1) + TGFbeta(1-j)
    end do
    TGFbeta(-b:0) = TGFbeta(L-b:L)
    TGFbeta(L+1:L+b+1) = TGFbeta(1:b+1)
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename, filename2
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    !    WRITE(filename2,'(A7,I5.5,A4)') './out/g', index, '.dat'
    open (unit = 11, file=filename, action="write")
    !    open (unit = 12, file=filename2, action="write")

    do i = 1, L+1
       do j = 1, H
          write(11, '(I5)', advance="no"), cmat(i,j)%type
       end do
       write(11, '(F15.2)', advance="no"), TGFbeta(i)
       write(11, *)
    end do
!!$    do i = 1, L
!!$       do j = 1, H
!!$          if (cmat(i,j)%type.eq.1) then
!!$             write(12, '(I10, (F15.5))'), i, cmat(i,j)%gene1
!!$          end if
!!$       end do
!!$    end do
    close(11)
    !    close(12)
  end subroutine output_to_file

  subroutine cell_event(i)
    use random
    implicit none
    integer, intent(in) :: i
    integer j, k, m
    real u, u1, p0
    real p_move(3)
    type(cell) new_cell
    call ran2(u)
    u = u*a(i)
    print *, 'event happen at ', i
    print *, 'u', u
    print *, 'a(i)', a(i)
    print *, 'npack(i)', npack(i-1:i+1)
    do j = 1, npack(i)
       u = u - v
       if ( u < 0 ) then
          if ( cmat(i,j)%type .eq. 3 ) then ! death
             !print *, 'die at', i, j, cmat(i,j)%type 
             !read(*,*)
             do k=j, H-1
                cmat(i, k) = cmat(i, k+1)
             end do
             TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) - D_TGFbeta
             npack(i) = npack(i) - 1
             return
          else if ( cmat(i,j)%type .eq. 1 ) then
             !print *, 'pro at', i, j, cmat(i,j)%type 
             !read(*,*)
             p0 = 0.2 + 0.6 / (1.0 + 0.01*TGFbeta(i))                
             call ran2(u1)
             if ( u1 < p0 ) then
                new_cell = cmat(i,j)
             else
                cmat(i,j)%type = 2
                new_cell = cmat(i,j) 
             end if
          else if ( cmat(i,j)%type .eq. 2 ) then
             !print *, 'pro at', i, j, cmat(i,j)%type 
             !read(*,*)
             call ran2(u1)
             if ( u1 < p1 ) then
                new_cell = cmat(i,j)
             else
                cmat(i,j)%type = 3
                new_cell = cmat(i,j)
                TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) + D_TGFbeta
             end if
          else
             print *, 'error 3'
             print *, npack(i), i, j
             print *, cmat(i,j)%type
             read(*,*)
          end if

          p_move(1) = 1.0/(npack(i-1)+1.0)
          p_move(2) = 0.0
          p_move(3) = 1.0/(npack(i+1)+1.0)
          call ran2(u1)
          u1 = u1*sum(p_move)
          do m = -1, 1
             u1 = u1 - p_move(m+2)
             if ( u1 < 0 ) then
                !print *, 'mov to', m, new_cell%type 
                !read(*,*)
                if ( cmat(i+m, j)%type .eq. 0 ) then
                   cmat(i+m, npack(i+m)+1) = new_cell                   
                else
                   do k=npack(i+m)+1, j+1, -1
                      cmat(i+m, k) = cmat(i+m, k-1)
                   end do
                   cmat(i+m, j) = new_cell
                   npack(i+m) = npack(i+m) + 1
                   if ( new_cell%type .eq. 3 ) then
                      TGFbeta(i+m-b:i+m+b) = TGFbeta(i+m-b:i+m+b) + D_TGFbeta
                   end if
                end if
                return
             end if
          end do
          return
       end if
    end do
    write(*,*) 'error 2'
    read(*,*)
  end subroutine cell_event

  subroutine cell_restack(i)
    use random
    implicit none
    integer, intent(in) :: i
    integer j
    type(cell) temp
    do j = npack(i), 2, -1
       if ( cmat(i,j)%type .eq. 1 ) then
          temp = cmat(i,j-1)
          cmat(i,j-1) = cmat(i,j)
          cmat(i,j) = temp
       end if
    end do
  end subroutine cell_restack

  subroutine cell_stat(t)
    implicit none

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
!!$  do i = -b, L+1+b
!!$     write(*, '(F8.2)', advance="no"), TGFbeta(i)
!!$  end do
!!$  write(*, *)
    write(100, '(5(F10.2))'), t, real(num_sc)/L, &
         real(num_tac)/L, real(num_tdc)/L, real(num_mc)/L
  end subroutine cell_stat

end module setting
