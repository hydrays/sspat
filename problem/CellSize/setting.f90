module setting
  integer :: NPool
  integer :: iseed
  real :: tpinc, tminc
  real :: tend, timestep
  real :: lambda1, gamma1, lambda2, gamma2, kappa
  real :: s0, newcell_delta 
  real :: age_critical, size_critical, omega, eta
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
	lambda1, gamma1, lambda2, gamma2, &
	s0, newcell_delta, kappa, size_critical, &
        age_critical, omega, eta

  type cell
     real Csize
     real Csize_old
     real Cage
     real mRNA
     real nRibsome
  end type cell

  type(cell), allocatable :: CellPool(:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'NPool = ', NPool
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'tminc = ', tminc
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'timestep = ', timestep
    write(*, '(a20, f10.2)'), 'lambda1 = ', lambda1
    write(*, '(a20, f10.2)'), 'gamma1 = ', gamma1
    write(*, '(a20, f10.2)'), 'lambda2 = ', lambda2
    write(*, '(a20, f10.2)'), 'gamma2 = ', gamma2
    write(*, '(a20, f10.2)'), 's0 = ', s0
    write(*, '(a20, f10.2)'), 'newcell_delta = ', newcell_delta
    write(*, '(a20, f10.2)'), 'kappa = ', kappa
    write(*, '(a20, f10.2)'), 'age_critical = ', age_critical
    write(*, '(a20, f10.2)'), 'size_critical = ', size_critical
    write(*, '(a20, f10.2)'), 'omega = ', omega
    write(*, '(a20, f10.2)'), 'eta = ', eta

    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'NPool,', NPool
    write(9, '(a20, i10)'), 'iseed,', iseed
    write(9, '(a20, f10.2)'), 'tpinc', tpinc
    write(9, '(a20, f10.2)'), 'tminc', tminc
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'timestep,', timestep
    write(9, '(a20, f10.2)'), 'lambda1,', lambda1
    write(9, '(a20, f10.2)'), 'gamma1,', gamma1
    write(9, '(a20, f10.2)'), 'lambda2,', lambda2
    write(9, '(a20, f10.2)'), 'gamma2,', gamma2
    write(9, '(a20, f10.2)'), 's0,', s0
    write(9, '(a20, f10.2)'), 'newcell_delta,', newcell_delta
    write(9, '(a20, f10.2)'), 'kappa,', kappa
    write(9, '(a20, f10.2)'), 'age_critical,', age_critical
    write(9, '(a20, f10.2)'), 'size_critical,', size_critical
    write(9, '(a20, f10.2)'), 'omega,', omega
    write(9, '(a20, f10.2)'), 'eta,', eta

    close(8)
    close(9)

  end subroutine read_xdata

  subroutine init_cell_pool()
    use random
    implicit none
    integer i

    allocate(CellPool(NPool))

    do i = 1, NPool
       CellPool(i)%Csize = s0
       CellPool(i)%Csize_old = s0
       CellPool(i)%Cage = 0.0
       CellPool(i)%mRNA = 0.0
       CellPool(i)%nRibsome = s0
    end do
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none
    integer, intent(in) :: index
    character(30) filename
    integer i

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, NPool
       write(11, '(6(F10.2))'), CellPool(i)%Csize, CellPool(i)%Csize_old, &
            CellPool(i)%Cage, CellPool(i)%mRNA, &
            CellPool(i)%nRibsome
    end do
    close(11)
  end subroutine output_to_file

  subroutine cell_division(i)
    use random
    implicit none
    integer, intent(in) :: i
    integer j
    real u, temp
    integer np
    type(cell) new_cell
    call normdev(0.0, newcell_delta, u)
    temp = CellPool(i)%Csize
    CellPool(i)%Csize = temp/2.0 + u
    new_cell%Csize = temp - CellPool(i)%Csize
    CellPool(i)%Csize_old = CellPool(i)%Csize
    new_cell%Csize_old = new_cell%Csize
    CellPool(i)%Cage = 0.0
    new_cell%Cage = 0.0
    CellPool(i)%mRNA = 0.0
    new_cell%mRNA = 0.0
    CellPool(i)%nRibsome = CellPool(i)%Csize
    new_cell%nRibsome = new_cell%Csize

    call ran2(u)
    np = ceiling(u*NPool)
    if (np .le. 0 .or. np > NPool) then
       print *, 'wrong np'
       read(*,*)
    end if
    CellPool(np) = new_cell
    
  end subroutine cell_division

  subroutine check_mitosis(size, age, m_flag)
    use random
    implicit none
    real, intent(in) :: size, age
    integer, intent(out) :: m_flag
    integer j
    real u, p
    m_flag = 0
    call ran2(u)
    p = omega*(size - size_critical)/size_critical + &
         (1-omega)*(age - age_critical)/age_critical
    p = eta * p * timestep
    if ( u < p ) then
       m_flag = 1
    end if
  end subroutine check_mitosis
  
  subroutine getrate(x, age, a)
    implicit none
    real, intent(in) :: x(lSpec)
    real, intent(in) :: age
    real, intent(out) :: a(lReac)

    a = 0.0
    a(1) = lambda1*kappa*age/(1.0+kappa*age)
    a(2) = gamma1*x(1)
    a(3) = lambda2*min(x(1), x(2))
    a(4) = gamma2*x(2)
  end subroutine getrate

  subroutine checkx(x, a, r, is_nag)
    implicit none
    real, intent(inout) :: x(lSpec)
    real, intent(inout) :: a(lReac)
    integer, intent(inout) :: r(lReac)
    integer,  intent(out) :: is_nag
    is_nag = 0
    if(any(x < 0.0) ) then
       is_nag = 1
       print *, x, a, r
       print *, 'nag!'
       read(*,*)
    end if
  end subroutine checkx
end module setting
