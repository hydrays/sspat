module setting
  integer Lbox
  integer :: iseed
  real :: tend, dt
  real :: alpha, beta
  integer :: pad
  
  namelist /xdata/ Lbox, tend, tpinc, iseed, dt, alpha, beta, pad

  integer, allocatable :: cmat(:,:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *) 'Control parameters...'
    write(*, '(a20, i10)') 'Lbox = ', Lbox
    write(*, '(a20, f10.2)') 'tpinc = ', tpinc
    write(*, '(a20, i10)') 'iseed = ', iseed
    write(*, '(a20, f10.2)') 'tend = ', tend
    write(*, '(a20, f10.2)') 'dt = ', dt
    write(*, '(a20, f10.2)') 'alpha = ', alpha
    write(*, '(a20, f10.2)') 'beta = ', beta
    write(*, '(a20, i10)') 'pad = ', pad

    open(9, file="out/control.csv")
    write(9, '(a20, a10)') 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)') 'Lbox,', Lbox
    write(9, '(a20, f10.210)') 'tpinc,', tpinc
    write(9, '(a20, i10)') 'iseed,', iseed
    write(9, '(a20, f10.2)') 'tend,', tend
    write(9, '(a20, f10.2)') 'dt,', dt
    write(9, '(a20, f10.2)') 'alpha,', alpha
    write(9, '(a20, f10.2)') 'beta,', beta
    write(9, '(a20, i10)') 'pad,', pad

    close(8)
    close(9)
  end subroutine read_xdata

  subroutine init()
    implicit none
    real u
    integer i, j
    !integer db_counter

    write(*, '(A)', advance='no') 'Initialize...'
    allocate(cmat(1:Lbox,1:Lbox))

    cmat = 0
    !db_counter = 0
    do i = pad+1, Lbox-pad
       do j = pad+1, Lbox-pad
          cmat(i,j) = 1
          !cmat(i,j) = db_counter
          !db_counter = db_counter + 1
       end do
    end do

    write(*, *) 'Done.'
  end subroutine init

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/c', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, Lbox
       do j = 1, Lbox-1
          write(11, '(I5, A2)', advance="no") cmat(i,j), ', '
       end do
       write(11, '(I5)') cmat(i,Lbox)
    end do
    close(11)
  end subroutine output_to_file

  subroutine padding()
    implicit none

    integer i, j
    integer L
    L = Lbox - 2*pad
    ! padding region 1
    do i = pad+1, Lbox-pad
       do j = 1, pad
          cmat(i,j) = cmat(i,j+L)
       end do
    end do
    ! padding region 2
    do i = Lbox-pad+1, Lbox
       do j = pad+1, Lbox-pad
          cmat(i,j) = cmat(i-L,j)
       end do
    end do
    ! padding region 3
    do i = pad+1, Lbox-pad
       do j = Lbox-pad+1, Lbox
          cmat(i,j) = cmat(i,j-L)
       end do
    end do
    ! padding region 4
    do i = 1, pad
       do j = pad+1, Lbox-pad
          cmat(i,j) = cmat(i+L,j)
       end do
    end do
    ! padding region 5
    do i = 1, pad
       do j = 1, pad
          cmat(i,j) = cmat(i+L,j+L)
       end do
    end do
    ! padding region 6
    do i = Lbox-pad+1, Lbox
       do j = 1, pad
          cmat(i,j) = cmat(i-L,j+L)
       end do
    end do
    ! padding region 7
    do i = Lbox-pad+1, Lbox
       do j = Lbox-pad+1, Lbox
          cmat(i,j) = cmat(i-L,j-L)
       end do
    end do
    ! padding region 8
    do i = 1, pad
       do j = Lbox-pad+1, Lbox
          cmat(i,j) = cmat(i+L,j-L)
       end do
    end do
  end subroutine padding
  
  subroutine get_energy(energy)
    implicit none
    real, intent(out) :: energy
    integer i, j
    energy = 0.0
    do i = pad+1, Lbox-pad
       do j = pad+1, Lbox-pad
          energy = energy - 0.5 * alpha * ( cmat(i, j) * cmat(i+1, j) &
               + cmat(i, j) * cmat(i-1, j) &
               + cmat(i, j) * cmat(i, j+1) &
               + cmat(i, j) * cmat(i, j-1) )
       end do
    end do
  end subroutine get_energy

  real function prob(E0, E1)
    implicit none
    real, intent(in) :: E0, E1

    prob = exp((E0 - E1) * beta)

  end function prob
  
end module setting
