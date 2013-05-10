module setting
  integer :: NPool
  integer :: iseed
  real :: tpinc
  real :: tend, timestep
  real :: lambda1, gamma1, lambda2, gamma2
  real :: alpha1, alpha2, kappa
  real :: s0, newcell_delta 
  real :: age_critical, size_critical, omega, eta
  integer, parameter :: lReac = 6
  integer, parameter :: lSpec = 3
  ! Reactions within one cell
  ! 1. DNA -> mRNA, a = lambda1*k*a/(1+k*a)
  ! 2. mRNA -> 0, a = gamma1*mRNA
  ! 3. mRNA + Ribsome -> ActiveRibsome, a = alpha1*mRNA*Ribsome 
  ! 4. ActiveRibsome -> mRNA + Ribsome, a = alpha2*ActiveRibsome
  ! 5. 0 -> Ribsome, a = lambda2*ActiveRibsome
  ! 6. Ribsome -> 0, a = gamma2*Ribsome
  integer, parameter, dimension(lSpec,lReac) :: nu = reshape( &
       ! 1   2   3   4   5
       source=(/ &
       (/01, 00, 00/), & !1
       (/-1, 00, 00/), & !2
       (/-1, 01, -1/), & !3
       (/01, -1, 01/), & !4
       (/00, 00, 01/), & !5
       (/00, 00, -1/) & !6
       /), shape = (/lSpec, lReac/) &
       )

  namelist /xdata/ NPool, iseed, tpinc,	tend, timestep, &
	lambda1, gamma1, lambda2, gamma2, &
	alpha1, alpha2, s0, newcell_delta, &
	kappa, size_critical, age_critical, omega, &
	eta

  type cell
     real Csize
     real Cage
     real Crate
     real mRNA
     real nActiveRibsome
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
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'timestep = ', timestep
    write(*, '(a20, f10.2)'), 'lambda1 = ', lambda1
    write(*, '(a20, f10.2)'), 'gamma1 = ', gamma1
    write(*, '(a20, f10.2)'), 'lambda2 = ', lambda2
    write(*, '(a20, f10.2)'), 'gamma2 = ', gamma2
    write(*, '(a20, f10.2)'), 'alpha1 = ', alpha1
    write(*, '(a20, f10.2)'), 'alpha2 = ', alpha2
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
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'timestep,', timestep
    write(9, '(a20, f10.2)'), 'lambda1,', lambda1
    write(9, '(a20, f10.2)'), 'gamma1,', gamma1
    write(9, '(a20, f10.2)'), 'lambda2,', lambda2
    write(9, '(a20, f10.2)'), 'gamma2,', gamma2
    write(9, '(a20, f10.2)'), 'alpha1,', alpha1
    write(9, '(a20, f10.2)'), 'alpha2,', alpha2
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
       CellPool(i)%Cage = 0.0
       CellPool(i)%Crate = 0.0
       CellPool(i)%mRNA = 0.0
       CellPool(i)%nRibsome = s0
       CellPool(i)%nActiveRibsome = 0.0
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
       write(11, '(6(F10.2))'), CellPool(i)%Csize, &
            CellPool(i)%Cage, CellPool(i)%Crate, CellPool(i)%mRNA, &
            CellPool(i)%nRibsome, CellPool(i)%nActiveRibsome
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
    CellPool(i)%Cage = 0.0
    new_cell%Cage = 0.0
    CellPool(i)%Crate = 0.0
    new_cell%Crate = 0.0
    CellPool(i)%mRNA = 0.0
    new_cell%mRNA = 0.0
    CellPool(i)%nActiveRibsome = 0.0
    new_cell%nActiveRibsome = 0.0
    CellPool(i)%nRibsome = CellPool(i)%Csize
    new_cell%Crate = new_cell%Csize

    call ran2(u)
    np = ceiling(u*NPool)
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
    a(3) = alpha1*x(1)*x(3)
    a(4) = alpha2*x(2)
    a(5) = lambda2*x(2)
    a(6) = gamma2*x(3)
  end subroutine getrate

  subroutine checkx(x, is_nag)
    implicit none
    real, intent(inout) :: x(lSpec)
    integer,  intent(out) :: is_nag
    is_nag = 0
    if(any(x < 0.0) ) then
       is_nag = 1
       print *, x
       print *, 'nag!'
       read(*,*)
    end if
  end subroutine checkx
end module setting
