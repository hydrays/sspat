module setting
  integer Lbox
  integer :: iseed
  real :: tend, dt
  real :: a0
  real :: beta, moblty
  integer :: R1
  integer :: nc
  real :: k0, phi0
  real :: lambda, difv, gamma
  
  namelist /xdata/ Lbox, tend, dt, a0, difv, &
       iseed, tpinc, R1, beta, nc, k0, &
       phi0, lambda, gamma, moblty

  type cell
     integer type
     ! type == 0 : empty splot
     ! type == 1 : M cell
     ! type == 2 : T cell
     ! type == 3 : activated T cell
  end type cell
  type(cell), allocatable :: cmat(:,:)
  real, allocatable ::  phi(:,:)
  real, allocatable ::  phi_old(:,:)
  real, allocatable :: p(:,:)
  real, allocatable :: a(:,:)
  
contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'Lbox = ', Lbox
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'dt = ', dt
    write(*, '(a20, f10.2)'), 'a0 = ', a0
    write(*, '(a20, f10.2)'), 'difv = ', difv
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, i10)'), 'R1 = ', R1
    write(*, '(a20, f10.2)'), 'beta = ', beta
    write(*, '(a20, f10.2)'), 'moblty = ', moblty
    write(*, '(a20, f10.2)'), 'k0 = ', k0
    write(*, '(a20, f10.2)'), 'phi0 = ', phi0
    write(*, '(a20, f10.2)'), 'lambda = ', lambda
    write(*, '(a20, f10.2)'), 'gamma = ', gamma
    write(*, '(a20, i10)'), 'nc = ', nc
    
    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'Lbox,', Lbox
    write(9, '(a20, i10)'), 'iseed,', iseed
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'dt,', dt
    write(9, '(a20, f10.2)'), 'a0,', a0
    write(9, '(a20, f10.2)'), 'difv,', difv
    write(9, '(a20, f10.2)'), 'tpinc,', tpinc
    write(9, '(a20, i10)'), 'R1,', R1
    write(9, '(a20, f10.2)'), 'beta,', beta
    write(9, '(a20, f10.2)'), 'moblty,', moblty
    write(9, '(a20, f10.2)'), 'k0,', k0
    write(9, '(a20, f10.2)'), 'phi0,', phi0
    write(9, '(a20, f10.2)'), 'lambda,', lambda
    write(9, '(a20, f10.2)'), 'gamma,', gamma
    write(9, '(a20, i10)'), 'nc,', nc
    
    close(8)
    close(9)
  end subroutine read_xdata
  
  subroutine init_cell_pool()
    implicit none
    real u
    integer i, j
    integer curb
    
    write(*, *), 'Initialize...'
    allocate(cmat(1:Lbox,1:Lbox))
    allocate(phi(1:Lbox, 1:Lbox))
    allocate(p(1:Lbox,1:Lbox))
    allocate(a(1:Lbox,1:Lbox))
    allocate(phi_old(1:Lbox,1:Lbox))

    do i = 1, Lbox
       do j = 1, Lbox
          cmat(i,j)%type = 0 ! No cell everywhere
          phi(i,j) = 1.0
          p(i,j) = 0.0
          a(i,j) = 0.0
       end do
    end do

    curb = 10
    ! randomly distribute M cells
    do i = curb, Lbox-curb, 4
       do j = curb, Lbox-curb, 4
          call random_number(u)
             cmat(i,j)%type = 1
       end do
    end do
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/c', index, '.dat'
    open (unit = 11, file=filename, action="write")
    WRITE(filename,'(A7,I5.5,A4)') './out/p', index, '.dat'
    open (unit = 21, file=filename, action="write")
    WRITE(filename,'(A9,I5.5,A4)') './out/phi', index, '.dat'
    open (unit = 31, file=filename, action="write")

    
    do i = 1, Lbox
       do j = 1, Lbox-1
          write(11, '(I5, A2)', advance="no"), cmat(i,j)%type, ', '
          write(21, '(F8.4, A2)', advance="no"), p(i,j), ', '
          write(31, '(F8.4, A2)', advance="no"), phi(i,j), ', '
       end do
       write(11, '(I5)'), cmat(i,Lbox)%type
       write(21, '(F8.4)'), p(i,j)
       write(31, '(F8.4)'), phi(i,j)
    end do
    close(11)
    close(21)
    close(31)
  end subroutine output_to_file

  subroutine cell_event(i, j)
    implicit none
    integer, intent(in) :: i, j
  end subroutine cell_event

  subroutine cell_stat()
    implicit none
  end subroutine cell_stat

  subroutine update_rate()
    implicit none
    integer i, j
    do i = 1, Lbox
       do j = 1, Lbox
          a(i,j) = a0*p(i, j)
       end do
    end do    
  end subroutine update_rate
  
  subroutine update_phi()
    implicit none
    integer i, j
    real A

    ! source
    do i = 1, Lbox
       do j = 1, Lbox
          if ( cmat(i,j)%type .eq. 3 ) then
             phi(i,j) = phi(i,j) + lambda*dt
          end if
       end do
    end do
    ! Neumann boundary condition (sort of)
    ! i = 1 or Lbox
    do j = 1, Lbox
       phi(1,j) = phi(2,j)
       phi(Lbox,j) = phi(Lbox-1,j)
    end do
    ! j = 1 or Lbox
    do i = 1, Lbox
       phi(i,1) = phi(i,2)
       phi(i,Lbox) = phi(i,Lbox-1)
    end do

    
    ! diffuse
    do i = 1, Lbox
       do j = 1, Lbox
          phi_old(i,j) = phi(i, j)
       end do
    end do
    do i = 2, Lbox-1
       do j = 2, Lbox-1
          phi(i, j) = phi_old(i,j) &
               + dt*difv*(phi_old(i-1,j)+phi_old(i+1,j)+phi_old(i,j-1)+phi_old(i,j+1)-4.0*phi_old(i,j))
       end do
    end do

    ! decay
    do i = 1, Lbox
       do j = 1, Lbox
          phi(i,j) = phi(i,j) - dt*gamma*phi(i,j)
       end do
    end do
    
    A = 0.0
    do i = 1, Lbox
       do j = 1, Lbox
          p(i, j) = 1.0/(1.0+exp(-k0*(phi(i, j) - phi0)))
          !p(i, j) = exp(k0*phi(i, j))
          A = A + p(i, j)
       end do
    end do

    A = A/(Lbox*Lbox)
    do i = 1, Lbox
       do j = 1, Lbox
          p(i, j) = p(i, j)/A
       end do
    end do
    
  end subroutine update_phi

end module setting
