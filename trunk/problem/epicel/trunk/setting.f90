module setting
  integer :: L
  integer :: brange, iseed
  real :: tend, p1, v, difv, mutv, tm
  real :: fdgain1, scstick, prelax, tpinc
  real :: bd10
  integer :: useomp, is64bit
  real :: timestep
  integer :: npar
  integer, parameter :: brange1 = 25

  namelist /xdata/ L, brange, tend, p1, v, difv, mutv, &
       fdgain1, scstick, prelax, iseed, tpinc, tm, HP0, HP1

  namelist /xdataomp/ useomp, is64bit, timestep, npar

  type cell
     integer type
     real HP
     real gene(3)
  end type cell
  type(cell), allocatable :: cmat(:,:)
  real, allocatable :: a(:)
  real, allocatable :: NT(:)
  real, allocatable :: NP(:)
  integer, allocatable :: npack(:)
  integer, allocatable :: TDC(:)
  integer, allocatable :: SC(:)
  integer, allocatable :: TAC(:)
  integer, allocatable :: MC(:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)
    read(8, nml=xdataomp)

    bd10 = 10.0/real(brange)

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'L = ', L
    write(*, '(a20, i10)'), 'brange = ', brange
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'p1 = ', p1
    write(*, '(a20, f10.2)'), 'v = ', v
    write(*, '(a20, f10.2)'), 'difv = ', difv
    write(*, '(a20, f10.2)'), 'mutv = ', mutv
    write(*, '(a20, f10.2)'), 'fdgain1 = ', fdgain1
    write(*, '(a20, f10.2)'), 'scstick = ', scstick
    write(*, '(a20, f10.2)'), 'prelax = ', prelax
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'tm = ', tm
    write(*, '(a20, f10.2)'), 'HP0 = ', HP0
    write(*, '(a20, f10.2)'), 'HP1 = ', HP1

    if (useomp.eq.1) then
       write(*, '(a)'), 'OpenMP parallel in use!'
       write(*, '(a20, I10)'), 'is64bit = ', is64bit
       write(*, '(a20, f10.2)'), 'timestep = ', timestep
       write(*, '(a20, I10)'), 'npar = ', npar
    end if

    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'L,', L
    write(9, '(a20, i10)'), 'brange,', brange
    write(9, '(a20, i10)'), 'iseed,', iseed
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'p1,', p1
    write(9, '(a20, f10.2)'), 'v,', v
    write(9, '(a20, f10.2)'), 'difv,', difv
    write(9, '(a20, f10.2)'), 'mutv,', mutv
    write(9, '(a20, f10.2)'), 'fdgain1,', fdgain1
    write(9, '(a20, f10.2)'), 'scstick,', scstick
    write(9, '(a20, f10.2)'), 'prelax,', prelax
    write(9, '(a20, f10.2)'), 'tpinc,', tpinc
    write(9, '(a20, f10.2)'), 'tm,', tm
    write(9, '(a20, i10)'), 'useomp,', useomp
    write(9, '(a20, i10)'), 'is64bit,', is64bit
    write(9, '(a20, f10.2)'), 'timestep,', timestep
    write(9, '(a20, I10)'), 'npar,', npar
    write(9, '(a20, f10.2)'), 'HP0,', HP0
    write(9, '(a20, f10.2)'), 'HP1,', HP1
    close(8)
    close(9)

  end subroutine read_xdata

  subroutine init_cell_pool()
    use random
    implicit none
    real u
    integer i, j

    allocate(cmat(-L:L,-L:L))
    do i = -L, L
       do j = -L, L
          cmat(i,j)%type = 0
       end do
    end do
    cmat(0, 0)%type = 1
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none
    integer, intent(in) :: index
    character(30) filename
    integer i, j
    
    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, L
       do j = 1, H
          write(11, '(I5)', advance="no"), cmat(i,j)%type
       end do
    end do
    close(11)
  end subroutine output_to_file

  subroutine cell_event(i, j)
    use random
    implicit none
    integer, intent(in) :: i, j
    integer k, m, shift_i, k1, k2
    real u, p0, TGFbeta
    real vr, vl
    type(cell) new_cell
    real u1, u2, Pa, Ta
    real pressure, Tpressure
    real a(3)
    
    ! a(1) --- diffusion
    ! a(2) --- division
    ! a(3) --- death
    a = 0
    if ( cmat(i,j)%type .eq. 0 ) then
       return
    else if ( cmat(i,j)%type .eq. 1 ) then
       a(1) = difv
       a(2) = a(1) + 1.0
       a(3) = a(2) + 0.0
       a = a * timestep
       call ran2(u)
       if ( u < a(1) ) then
          ! diffusion
          
       else if ( u < a(2) ) then
          ! division
       else if (u < a(3) ) then
          ! death
          print *, 'should not happen'
       else
          ! do nothing
       end if
    end if

  end subroutine cell_event

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

  subroutine Perodic_BC(k)
    implicit none
    integer, intent(in) :: k

    if ( k .eq. 1 ) then
       cmat(L, :) = cmat(0, :)
       cmat(L+1, :) = cmat(1, :)
       npack(L) = npack(0)
       npack(L+1) = npack(1)
       SC(L) = SC(0)
       SC(L+1) = SC(1)
       TAC(L) = TAC(0)
       TAC(L+1) = TAC(1)
       TDC(L) = TDC(0)
       TDC(L+1) = TDC(1)
       MC(L) = MC(0)
       MC(L+1) = MC(1)
    end if
    if ( k .eq. 2 ) then
       cmat(L+1, :) = cmat(1, :)
       npack(L+1) = npack(1)
       TDC(L+1) = TDC(1)
       SC(L+1) = SC(1)
       TAC(L+1) = TAC(1)
       MC(L+1) = MC(1)
    end if
    if ( k .eq. L ) then
       cmat(1, :) = cmat(L+1, :)
       cmat(0, :) = cmat(L, :)
       npack(1) = npack(L+1)
       npack(0) = npack(L)
       SC(1) = SC(L+1)
       SC(0) = SC(L)
       TAC(1) = TAC(L+1)
       TAC(0) = TAC(L)
       TDC(1) = TDC(L+1)
       TDC(0) = TDC(L)
       MC(1) = MC(L+1)
       MC(0) = MC(L)
    end if
    if ( k .eq. L-1 ) then
       cmat(0, :) = cmat(L, :)
       npack(0) = npack(L)
       SC(0) = SC(L)
       TAC(0) = TAC(L)
       TDC(0) = TDC(L)
       MC(0) = MC(L)
    end if
  end subroutine Perodic_BC

  subroutine Next_Reaction(k, tau)
    implicit none
    integer, intent(out) :: k
    real, intent(out) :: tau
    real tau_temp
    integer i

    tau = huge(0.0)
    k = 0
    do i = 1, L
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
    if ( tau < 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, 'NP-NT', NP - NT
       read(*,*)
    end if
  end subroutine Next_Reaction

  subroutine Update_Rate(i)
    implicit none
    integer, intent(in) :: i
    real vr, vl
    vr = max(0.0, difv*real(npack(i)-npack(i+1)))
    vl = max(0.0, difv*real(npack(i)-npack(i-1)))
    !vr = difv*real(npack(i))
    !vl = difv*real(npack(i))
    if (npack(i).eq.0) then
       !       print *, vr, vl
    end if
    !a(i) = vr + vl + (SC(i)+TAC(i)+TDC(i)+MC(i))*(v+mutv)
    a(i) = vr + vl + npack(i)*(v+mutv)
  end subroutine Update_Rate
end module setting
