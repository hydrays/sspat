module setting
  integer, parameter :: L = 400
  integer, parameter :: H = 200
  real, parameter :: b = 8.0
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 10000.0
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
    type(cell) new_cell
    call ran2(u)
    u = u*a(i)
    if (npack(i).ne.0) then
       vr = max(0.0, 200.0*real(npack(i)-npack(i+1)))
       vl = max(0.0, 200.0*real(npack(i)-npack(i-1)))
    else
       vr = 0.0
       vl = 0.0
    end if
!!$    print *, 'event happen at ', i
!!$    print *, 'u', u
!!$    print *, 'a(i)', a(i)
!!$    print *, 'npack(i)', npack(i-1:i+1)
!!$    print *, 'vr, vl', vr, vl
    do j = 1, npack(i)
       u = u - v
       if ( u < 0 ) then
          if ( cmat(i,j)%type .eq. 1 ) then
             p0 = 0.2 + 0.6 / (1.0 + 0.01*TGFbeta(i))
             ! division
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             call ran2(u1)
             if ( u1 < p0 ) then
                ! SC -> 2SC
                cmat(i,j+1) = cmat(i,j)
             else
                ! SC -> 2TAC
                cmat(i,j)%type = 2
                cmat(i,j+1) = cmat(i,j) 
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 2 ) then
             ! division
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             call ran2(u1)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i, j+1) = cmat(i,j)
             else
                ! TAC -> 2TDC
                cmat(i, j)%type = 3
                cmat(i, j+1) = cmat(i,j)
                TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) + 2.0*D_TGFbeta
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 3 ) then
             ! death
             do k=j, H-1
                cmat(i, k) = cmat(i, k+1)
             end do
             TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) - D_TGFbeta
             npack(i) = npack(i) - 1
          else
             ! do nothing
          end if
          return
       end if
    end do
    u = u - vr
    if ( u < 0 ) then
       call ran2(u1)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 4', j, u1
          read(*,*)
       end if
       !print *, 'move right at height j', i, j

       if ( cmat(i, j)%type .eq. 1 ) then
          call ran2(u1)
          if (u1 < 0.98) then
             return
          end if
       end if

       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i + 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1
       if (new_cell%type .eq. 3) then
          TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) - D_TGFbeta
          TGFbeta(m-b:m+b) = TGFbeta(m-b:m+b) + D_TGFbeta
       end if
       return
    end if
    u = u - vl
    if ( u < 0 ) then
       call ran2(u1)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 5', j, u1, npack(i)
          read(*,*)
       end if
       !print *, 'move left at height j', i, j
       if ( cmat(i, j)%type .eq. 1 ) then
          call ran2(u1)
          if (u1 < 0.98) then
             return
          end if
       end if
       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i - 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1
       if (new_cell%type .eq. 3) then
          TGFbeta(i-b:i+b) = TGFbeta(i-b:i+b) - D_TGFbeta
          TGFbeta(m-b:m+b) = TGFbeta(m-b:m+b) + D_TGFbeta
       end if
       return
    else
!!$       print *, 'u', u
!!$       print *, 'vl', vl
       write(*,*) 'error 2'
       read(*,*)
    end if
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
