module setting
  integer :: NPool
  integer :: iseed
  real :: tpinc, tminc
  real :: tend, timestep
  real :: rho
  real :: lambda1, gamma1, lambda2, gamma2, kappa
  real :: s0, newcell_delta 
  real :: mp1, mp2, mp3
  real :: tc
  integer :: NCollect, NCollect2
  integer :: NSize

  real, parameter :: ctime = 14.0
  integer, parameter :: lReac = 4
  integer, parameter :: lSpec = 2
  ! Reactions within one cell
  ! 1. DNA -> mRNA, a = lambda1*k*a/(1+k*a)
  ! 2. mRNA -> 0, a = gamma1*mRNA
  ! 3. 0 -> Ribsome, a = lambda2*ActiveRibsome
  ! 4. Ribsome -> 0, a = gamma2*Ribsome
  integer, parameter, dimension(lSpec,lReac) :: nu = reshape( &
       ! 1   2
       source=(/ &
       (/01, 00/), & !1
       (/-1, 00/), & !2
       (/00, 01/), & !3
       (/00, -1/) & !4
       /), shape = (/lSpec, lReac/) &
       )

  namelist /xdata/ NPool, iseed, tpinc,	tminc, tend, timestep, &
       rho, lambda1, gamma1, lambda2, gamma2, &
       s0, newcell_delta, kappa, NCollect, NCollect2, &
       tc, NSize

  namelist /MitosisPara/ mp1, mp2!, mp3

  type cell
     ! size, age, mRNA, event
     real State(4)
     ! cellsize_old1, cellsize_old2
     real AuxVariable(2)
     ! growth rate
     real Pheno(8)
     real mitmet
  end type cell

  type(cell), allocatable :: CellPool(:)
  type(cell), allocatable :: MitosisPool(:)
  type(cell), allocatable :: NewbornPool(:)
  type(cell), allocatable :: NewbornPool2(:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)
    read(8, nml=MitosisPara)

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'NPool = ', NPool
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'tminc = ', tminc
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'timestep = ', timestep
    write(*, '(a20, f10.2)'), 'rho = ', rho
    write(*, '(a20, f10.2)'), 'lambda1 = ', lambda1
    write(*, '(a20, f10.2)'), 'gamma1 = ', gamma1
    write(*, '(a20, f10.2)'), 'lambda2 = ', lambda2
    write(*, '(a20, f10.2)'), 'gamma2 = ', gamma2
    write(*, '(a20, f10.2)'), 's0 = ', s0
    write(*, '(a20, f10.2)'), 'newcell_delta = ', newcell_delta
    write(*, '(a20, f10.2)'), 'kappa = ', kappa
    write(*, '(a20, I10)'), 'NCollect = ', NCollect
    write(*, '(a20, I10)'), 'NCollect2 = ', NCollect2
    write(*, '(a20, f10.2)'), 'tc ', tc
    write(*, '(a20, I10)'), 'NSize ', NSize

    write(*, *)

    write(*, *), 'Mitosis Parameters'
    write(*, '(a20, f10.2)'), 'mp1', mp1
    write(*, '(a20, f10.2)'), 'mp2', mp2
!    write(*, '(a20, f10.2)'), 'mp3', mp3

    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'NPool,', NPool
    write(9, '(a20, i10)'), 'iseed,', iseed
    write(9, '(a20, f10.2)'), 'tpinc,', tpinc
    write(9, '(a20, f10.2)'), 'tminc,', tminc
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'timestep,', timestep
    write(9, '(a20, f10.2)'), 'rho,', rho
    write(9, '(a20, f10.2)'), 'lambda1,', lambda1
    write(9, '(a20, f10.2)'), 'gamma1,', gamma1
    write(9, '(a20, f10.2)'), 'lambda2,', lambda2
    write(9, '(a20, f10.2)'), 'gamma2,', gamma2
    write(9, '(a20, f10.2)'), 's0,', s0
    write(9, '(a20, f10.2)'), 'newcell_delta,', newcell_delta
    write(9, '(a20, f10.2)'), 'kappa,', kappa
    write(9, '(a20, I10)'), 'NCollect,', NCollect
    write(9, '(a20, I10)'), 'NCollect2,', NCollect2
    write(9, '(a20, f10.2)'), 'tc,', tc
    write(9, '(a20, I10)'), 'NSize,', NSize
    write(9, '(a20, f10.2)'), 'mp1,', mp1
    write(9, '(a20, f10.2)'), 'mp2,', mp2
!    write(9, '(a20, f10.2)'), 'mp3,', mp3

    close(8)
    close(9)

    print *, "Multiply lambda1 by rho"
    lambda1 = lambda1*rho
  end subroutine read_xdata

  subroutine init_cell_pool()
    use random
    implicit none
    integer i
    real u

    allocate(CellPool(NPool))

    allocate(MitosisPool(NPool))
    allocate(NewbornPool(NPool))
    allocate(NewbornPool2(NPool))

    do i = 1, NPool
       call ran2(u)
       CellPool(i)%Pheno(1) = lambda1
       CellPool(i)%Pheno(2) = gamma1
       CellPool(i)%Pheno(3) = lambda2 + 0.5*lambda2*(2.0*u - 1.0)
       ! if ( u < 0.5 ) then
       !    CellPool(i)%Pheno(3) = 0.12
       ! else
       !    CellPool(i)%Pheno(3) = 0.08
       ! end if       
       !CellPool(i)%Pheno(4) = ctime - log(2.0)/CellPool(i)%Pheno(3)
       !print *, CellPool(i)%mitmet

       s0 = 500.0/(0.35+CellPool(i)%Pheno(3))
       CellPool(i)%State(1) = s0
       CellPool(i)%AuxVariable(1) = s0
       CellPool(i)%AuxVariable(2) = s0
       CellPool(i)%State(2) = 0.0
       CellPool(i)%State(3) = s0
       call ran2(u)
       CellPool(i)%mitmet = log(u)
       CellPool(i)%State(4) = -1.0       
       
       CellPool(i)%Pheno(4) = max(0.0, ctime - log(2.0)/CellPool(i)%Pheno(3))
       !CellPool(i)%Pheno(4) = max(0.0, ctime - (1.0/CellPool(i)%Pheno(3))*&
       !     log(1000.0/(CellPool(i)%State(1)*(CellPool(i)%Pheno(3)+0.35))))

    end do
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none
    integer, intent(in) :: index
    character(30) filename
    integer i

    if ( index .ge. 0 ) then
       WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
       open (unit = 11, file=filename, action="write")

       do i = 1, NPool
          write(11, '(6(F16.2))'), &
               CellPool(i)%State(1), &
               CellPool(i)%AuxVariable(1), &
               CellPool(i)%AuxVariable(2), &
               CellPool(i)%State(2), &
               CellPool(i)%Pheno(3), &
               CellPool(i)%State(4)
       end do
       close(11)
    end if
  end subroutine output_to_file

  subroutine output_to_file_newborn(index)
    implicit none
    integer, intent(in) :: index
    character(30) filename
    integer i

    WRITE(filename,'(A7,I5.5,A4)') './out/n', index, '.dat'
    open (unit = 11, file=filename, action="write")
    do i = 1, NPool
       write(11, '((F16.2))'), NewbornPool(i)%State(1)
    end do
    close(11)
  end subroutine output_to_file_newborn

  subroutine cell_division(i)
    use random
    implicit none
    integer, intent(in) :: i
    integer j
    real u, temp
    integer np
    type(cell) new_cell
    call normdev(0.0, newcell_delta, u)
    !call normdev(0.0, newcell_delta*CellPool(i)%State(1)/500, u)
    !write(17, '(6(F16.2))'), u
    u = min(u, 0.9*CellPool(i)%State(1))
    u = max(u, -0.9*CellPool(i)%State(1))
    temp = CellPool(i)%State(1)
    CellPool(i)%State(1) = (temp + u)/2.0
    new_cell%State(1) = temp - CellPool(i)%State(1)
    if ( CellPool(i)%State(1) * new_cell%State(1) .le. 0.0 ) then
       print *, 'Cell born to be negative size, abort ...'
       print *, newcell_delta, u, CellPool(i)%State(1)
    end if
    CellPool(i)%AuxVariable(1) = CellPool(i)%State(1)
    new_cell%AuxVariable(1) = new_cell%State(1)
    CellPool(i)%AuxVariable(2) = CellPool(i)%State(1)
    new_cell%AuxVariable(2) = new_cell%State(1)
    CellPool(i)%State(2) = 0.0
    new_cell%State(2) = 0.0
    new_cell%State(3) = 0.0*CellPool(i)%State(3)*(1.0-CellPool(i)%State(1)/temp)
    CellPool(i)%State(3) = 0.0*CellPool(i)%State(3)*CellPool(i)%State(1)/temp
    call ran2(u)
    CellPool(i)%mitmet = log(u)
    call ran2(u)
    new_cell%mitmet = log(u)
    CellPool(i)%State(4) = -1.0
    new_cell%State(4) = -1.0

    new_cell%Pheno = CellPool(i)%Pheno
    !CellPool(i)%Pheno(4) = max(0.0, ctime - (1.0/CellPool(i)%Pheno(3))*&
    !     log(1000.0/(CellPool(i)%State(1)*(CellPool(i)%Pheno(3)+0.35))))
    !new_cell%Pheno(4) = max(0.0, ctime - (1.0/new_cell%Pheno(3))*&
    !     log(1000.0/(new_cell%State(1)*(new_cell%Pheno(3)+0.35))))

    call ran2(u)
    np = ceiling(u*NPool)
    if (np .le. 0 .or. np > NPool) then
       print *, 'wrong np'
       !read(*,*)
    end if
    CellPool(np) = new_cell

  end subroutine cell_division

  ! subroutine cell_division_syn(i)
  !   use random
  !   implicit none
  !   integer, intent(in) :: i
  !   integer j
  !   real u, temp
  !   integer np
  !   type(cell) new_cell
  !   call normdev(0.0, newcell_delta, u)
  !   !call normdev(0.0, newcell_delta*NewbornPool(i)%State(1)/500, u)
  !   u = min(u, 0.9*NewbornPool(i)%State(1))
  !   u = max(u, -0.9*NewbornPool(i)%State(1))
  !   temp = NewbornPool(i)%State(1)
  !   NewbornPool(i)%State(1) = (temp + u)/2.0
  !   new_cell%State(1) = temp - NewbornPool(i)%State(1)
  !   if ( NewbornPool(i)%State(1) * new_cell%State(1) .le. 0.0 ) then
  !      print *, 'Cell born to be negative size, abort ...'
  !   end if
  !   NewbornPool(i)%AuxVariable(1) = NewbornPool(i)%State(1)
  !   new_cell%AuxVariable(1) = new_cell%State(1)
  !   NewbornPool(i)%AuxVariable(2) = NewbornPool(i)%State(1)
  !   new_cell%AuxVariable(2) = new_cell%State(1)
  !   NewbornPool(i)%State(2) = 0.0
  !   new_cell%State(2) = 0.0
  !   NewbornPool(i)%State(3) = 0.0
  !   new_cell%State(3) = 0.0
  !   NewbornPool(i)%nRibsome = rho*NewbornPool(i)%State(1)
  !   new_cell%nRibsome = rho*new_cell%State(1)

  !   call ran2(u)
  !   np = ceiling(u*NCollect)
  !   if (np .le. 0 .or. np > NCollect) then
  !      print *, 'wrong np'
  !      !read(*,*)
  !   end if
  !   NewbornPool(np) = new_cell
  ! end subroutine cell_division_syn

  subroutine check_mitosis(i, m_flag)
    use random
    implicit none
    integer, intent(in) :: i
    integer, intent(out) :: m_flag
    real event, size, age
    real p, p1, p2
    real lambda

    ! Scheme event
    ! event = CellPool(i)%State(4)
    ! p = mp2*(1.0+tanh(10.0*(event/(mp1*8000.0) - 1.0)))/2.0

    ! ! Scheme event 2
    ! event = CellPool(i)%State(4)
    ! if ( event > mp1*8000.0 ) then
    !    p = mp2
    ! else
    !    p = 0.0
    ! end if

    ! ! Scheme size + time
    ! size = CellPool(i)%State(1)
    ! age = CellPool(i)%State(2)
    ! if ( age > mp3*8.0 ) then
    !    !p1 = 0.2*max(0.0, (age/6.5 - 1.0))
    !    p1 = mp2
    ! else
    !    p1 = 0.0
    ! end if
    ! if (size > mp1*1200.0 ) then
    !    !p = mp2*(1.0+tanh(2.0*(size/(mp1*1200.0) - 1.0)))/2.0
    !    !p = mp2*(max(0.0, (size/(mp1*1200.0) - 1.0)))**3.0
    !    p2 = mp2
    ! else
    !    p2 = 0.0
    ! end if
    ! p = p1 + p2

    ! ! Scheme time only
    ! age = CellPool(i)%State(2)
    ! if ( age > mp1 ) then
    !     p = mp2*max(0.0, (age/(mp1*8.0) - 1.0))
    !  else
    !     p = 0.0
    ! end if

    ! Scheme time only
    size = CellPool(i)%State(1)
    lambda = 0.35+CellPool(i)%Pheno(3)
    if ( size*lambda > 1000.0 ) then
        p = 5.0
     else
        p = 0.0
    end if

    CellPool(i)%mitmet = CellPool(i)%mitmet + p*timestep
    if ( CellPool(i)%mitmet .ge. 0.0 ) then
       m_flag = 1
       !print *, 'cell divide at', size, age
       !read(*,*)
    else
       m_flag = 0
    end if
  end subroutine check_mitosis

  subroutine getrate(i, x, age, a)
    implicit none
    real, intent(in) :: x(lSpec)
    real, intent(in) :: age
    real, intent(out) :: a(lReac)
    integer, intent(in) :: i
    a(1) = lambda1*((kappa*age)**4)/(1.0+((kappa*age)**4))
    a(2) = gamma1*x(1)
    !a(3) = lambda2*min(x(1), x(2))
    !a(3) = lambda2*x(2)
    a(3) = CellPool(i)%Pheno(3)*x(1)
    a(4) = gamma2*x(2)
  end subroutine getrate

  subroutine getderivative(i, t, x, yp)
    implicit none
    real, intent(in) :: t
    real, intent(in) :: x(lSpec)
    real, intent(out) :: yp(lSpec)
    integer, intent(in) :: i
    yp(1) = lambda1*((kappa*t)**4)/(1.0+((kappa*t)**4)) - gamma1*x(1)
    !yp(2) = max(0.0, lambda2*min(x(1), x(2)) - gamma2*x(2))
    if ( CellPool(i)%State(2) > CellPool(i)%Pheno(4) ) then
       yp(2) = CellPool(i)%Pheno(3)*x(1)
    else
       yp(2) = 0.0
    end if
  end subroutine getderivative

  subroutine checkx(x, a, r, is_nag)
    implicit none
    real, intent(inout) :: x(lSpec)
    real, intent(inout) :: a(lReac)
    integer, intent(inout) :: r(lReac)
    integer,  intent(out) :: is_nag
    is_nag = 0
    if(any(x < 0.0) ) then
       is_nag = 1
       print *, x
       print *, a
       print *, r
       print *, 'nag!'
       !read(*,*)
    end if
  end subroutine checkx

end module setting
