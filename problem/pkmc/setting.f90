module setting
  integer, parameter :: L = 2000
  integer, parameter :: H = 200
  real, parameter :: b = 8.0
  real, parameter :: tend = 10000.0
  real, parameter :: p1 = 0.3
  real, parameter :: v = 1.0
  real, parameter :: stem_fix = 0.96
  type cell
     integer type
     real gene1
  end type cell
  type(cell) cmat(0:L+1,H)
  real a(1:L)
  real NT(1:L)
  real NP(1:L)
  integer npack(0:L+1)
  integer TDC(0:L+1)
  integer SC(0:L+1)
  integer TAC(0:L+1)

contains
  subroutine init_cell_pool()
    use par_zig_mod
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
    TDC = 0
    do i = 1, L
       do j = 1, H
          if ( cmat(i,j)%type .ne. 0 ) then
             npack(i) = npack(i) + 1
          end if
          if ( cmat(i,j)%type .eq. 3 ) then
             TDC(i) = TDC(i) + 1
          end if
       end do
    end do
    npack(0) = npack(L)
    npack(L+1) = npack(1)
    TDC(0) = TDC(L)
    TDC(L+1) = TDC(1)

    NP = 0.0
    NT = 0.0
    do i = 1, L
       call random_number(u)
       u = -log(u)
       call Update_Rate(i)
       NP(i) = u
    end do

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
       write(11, '(I6)', advance="no"), TDC(i)
       write(11, '(I6)', advance="no"), TAC(i)
       write(11, '(I6)', advance="no"), SC(i)
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

  subroutine cell_event(i, kpar)
    use par_zig_mod
    implicit none
    integer, intent(in) :: i, kpar
    integer j, k, m, shift_i
    real u, u1, p0, TGFbeta
    real vr, vl
    type(cell) new_cell
    u = par_uni(kpar)
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
             TGFbeta = 0.0
             do k = -b, b
                shift_i = k + i
                if ( shift_i .le. 0 ) then
                   shift_i = shift_i + L
                else if ( shift_i > L ) then
                   shift_i = shift_i - L
                end if
                TGFbeta = TGFbeta + TDC(shift_i)*exp(-real(abs(k))/b)
             end do
             p0 = 0.2 + 0.6 / (1.0 + 0.01*TGFbeta)
             !print *, 'p0', p0
             ! division
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             u1 = par_uni(kpar)
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
             u1 = par_uni(kpar)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i, j+1) = cmat(i,j)
             else
                ! TAC -> 2TDC
                cmat(i, j)%type = 3
                cmat(i, j+1) = cmat(i,j)
                TDC(i) = TDC(i) + 2
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 3 ) then
             ! death
             do k=j, H-1
                cmat(i, k) = cmat(i, k+1)
             end do
             TDC(i) = TDC(i) - 1
             npack(i) = npack(i) - 1
          else
             ! do nothing
          end if
          return
       end if
    end do
    u = u - vr
    if ( u < 0 ) then
       u1 = par_uni(kpar)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 4', j, u1
          read(*,*)
       end if

       if ( cmat(i, j)%type .eq. 1 ) then
          u1 = par_uni(kpar)
          if (u1 < stem_fix) then
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
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       end if
       return
    end if
    u = u - vl
    if ( u < 0 ) then
       u1 = par_uni(kpar)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 5', j, u1, npack(i)
          read(*,*)
       end if
       !print *, 'move left at height j', i, j
       if ( cmat(i, j)%type .eq. 1 ) then
          u1 = par_uni(kpar)
          if (u1 < stem_fix) then
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
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       end if
       return
    else
       print *, 'event at', i
       print *, 'npack', npack(i-1:i+1)
       print *, 'cell', cmat(i, 1:npack(i))%type
       print *, 'u', u
       print *, 'vl', vl
       print *, 'a', a(i)
       write(*,*) 'error 2'
       read(*,*)
    end if
  end subroutine cell_event

  subroutine cell_restack(i)
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
    integer inum_sc, inum_tac, inum_tdc, inum_mc

    num_sc = 0
    num_tac = 0
    num_tdc = 0
    num_mc = 0
    do i = 1, L
       inum_sc = 0
       inum_tac = 0
       inum_tdc = 0
       inum_mc = 0
       do j = 1, H
          if (cmat(i,j)%type.eq.1) then
             num_sc = num_sc + 1
             inum_sc = inum_sc + 1
          elseif (cmat(i,j)%type.eq.2) then
             num_tac = num_tac + 1
             inum_tac = inum_tac + 1
          elseif (cmat(i,j)%type.eq.3) then
             num_tdc = num_tdc + 1
             inum_tdc = inum_tdc + 1
          elseif (cmat(i,j)%type.eq.4) then
             num_mc = num_mc + 1
             inum_mc = inum_mc + 1
          end if
       end do
       SC(i) = inum_sc
       TAC(i) = inum_tac
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

  subroutine Perodic_BC(k)
    implicit none
    integer, intent(in) :: k

    if ( k .eq. 1 ) then
       cmat(L, :) = cmat(0, :)
       cmat(L+1, :) = cmat(1, :)
       npack(L) = npack(0)
       npack(L+1) = npack(1)
       TDC(L) = TDC(0)
       TDC(L+1) = TDC(1)
    end if
    if ( k .eq. 2 ) then
       cmat(L+1, :) = cmat(1, :)
       npack(L+1) = npack(1)
       TDC(L+1) = TDC(1)
    end if
    if ( k .eq. L ) then
       cmat(1, :) = cmat(L+1, :)
       cmat(0, :) = cmat(L, :)
       npack(1) = npack(L+1)
       npack(0) = npack(L)
       TDC(1) = TDC(L+1)
       TDC(0) = TDC(L)
    end if
    if ( k .eq. L-1 ) then
       cmat(0, :) = cmat(L, :)
       npack(0) = npack(L)
       TDC(0) = TDC(L)
    end if
  end subroutine Perodic_BC

  subroutine Next_Reaction(k, tau, ilow, iup)
    implicit none
    integer, intent(out) :: k
    integer, intent(in) :: ilow, iup
    real, intent(out) :: tau
    real tau_temp
    integer i

    tau = huge(0.0)
    k = 0
    do i = ilow, iup
       if ( a(i) > 0.0 ) then
          tau_temp = ( NP(i) - NT(i) ) / a(i)
          if ( tau_temp < tau) then
             tau = tau_temp
             k = i
          end if
       else
!          print *, 'no reaction -------------- ', i, k, a(i)
!          read(*,*)
       end if
    end do
!    print *, 'k', k
!    print *, 'a', a
    if ( k .le. 0 ) then
       write(*,*), 'error'
       read(*,*)
    end if
    if ( tau .le. 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, 'NP-NT', NP - NT
       read(*,*)
    end if
  end subroutine Next_Reaction

  subroutine Update_Rate(i)
    implicit none
    integer, intent(in) :: i
    real vr, vl
    vr = max(0.0, 200.0*real(npack(i)-npack(i+1)))
    vl = max(0.0, 200.0*real(npack(i)-npack(i-1)))
    if (npack(i).eq.0) then
!       print *, vr, vl
    end if
    a(i) = vr + vl + npack(i)*v
  end subroutine Update_Rate

end module setting
