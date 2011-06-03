module paras  
  use random
  implicit none
 
  real(kind=8), parameter :: cell_life_length = 1.0
  real(kind=8), parameter :: lambda0 = 10.0
  integer num, num_cell_can_born, next_die_index, next_born_index
  integer num_of_tdc, num_old
  real(kind=8) next_die_time, next_born_time
  integer next_die_index_old, next_born_index_old
  real(kind=8) next_die_time_old, next_born_time_old
  real(kind=8)  lambda, p, t, p0, p1
  real(kind=8), parameter :: tend=500
  real(kind=8), parameter :: L1=2000
!  real(kind=8), parameter :: L2=1000
  real(kind=8), parameter :: k=0.1
  
  type cell
     integer celltype
     real(kind=8) cellclock
  end type cell

contains
  
  subroutine output_to_file(cell_pool, num, index)
    implicit none
    type (cell), intent(in) :: cell_pool(:)
    integer, intent(in) :: index, num
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, num
       write(11, '(I2)', advance="no"), cell_pool(i)%celltype
    end do
    write(11, '(3(F6.2))', advance="no"), lambda, p, t
    if (next_die_time < next_born_time) then 
       write(11, '(I8)', advance="no"), next_die_index
       write(11, '(I3)', advance="no"), -1
    else
       write(11, '(I8)', advance="no"), next_born_index
       write(11, '(I3)', advance="no"), 1
    end if
    close(11)
  end subroutine output_to_file

  subroutine cell_count(cell_pool, num)
    implicit none
    type (cell), intent(in) :: cell_pool(:)
    integer, intent(in) :: num
    integer i
    integer num_sc, num_tac, num_tdc

    num_sc = 0
    num_tac = 0
    num_tdc = 0
    do i = 1, num
       if ( cell_pool(i)%celltype .eq. 1 ) then 
          num_sc = num_sc + 1
       elseif ( cell_pool(i)%celltype .eq. 2 ) then 
          num_tac = num_tac + 1
       elseif ( cell_pool(i)%celltype .eq. 3 ) then 
          num_tdc = num_tdc + 1
       else
          print *, 'woring'
          read(*,*)
       end if
    end do
    write(*, '(F10.2, 4I6, F10.4)'), t, num_sc, num_tac, num_tdc, num, p0
  end subroutine cell_count

  subroutine pool_print(unit, cell_pool, num)
    implicit none
    integer, intent(in) :: unit
    type (cell), intent(in) :: cell_pool(:)
    integer, intent(in) :: num
    
    integer i
    do i = 1, num
       write(unit, *), cell_pool(i)%celltype
    end do
  end subroutine pool_print

  subroutine pool_init(cell_pool, num)
    implicit none
    type (cell), intent(out) :: cell_pool(:)
    integer, intent(out) :: num
    real(kind=8) u
    real(kind=8) r

    integer i
!!$    num = 20
!!$    cell_pool(1:5)%celltype = 1
!!$    cell_pool(6:10)%celltype = 2
!!$    cell_pool(11:20)%celltype = 3
!!$    do i = 1, num
!!$       call ran2(u)
!!$       cell_pool(i)%cellclock = u*cell_life_length
!!$    end do

    num = 10
    cell_pool(1:10)%celltype = 1
    cell_pool(1:10)%cellclock = cell_life_length    
  end subroutine pool_init
  
end module paras
