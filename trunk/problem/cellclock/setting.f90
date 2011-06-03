module paras  
 
  real(kind=8), parameter :: lifespan(3) = (/1.0, 1.0, 1.0/)
  real(kind=8), parameter :: lambda0 = 10.0
  integer num, num_cell_can_born, die_index, born_index
  integer num_of_tdc, num_old
  real(kind=8) die_time, born_time
  !integer die_index_old, born_index_old
  !real(kind=8) die_time_old, born_time_old
  real(kind=8)  p(3), t
  real(kind=8), parameter :: tend=500
  real(kind=8), parameter :: L1=2000
  real(kind=8), parameter :: k=0.01
  
  type cell
     integer type
     real(kind=8) clock
  end type cell

contains

  subroutine getp
    implicit none    
    p(1) = 0.2 + (0.8 - 0.2)/(1.01 + k*real(num_of_tdc)/real(L1))
    !p(1) = 1.0/(1.01 + k*real(num_of_tdc)/real(L1))
    p(2) = 0.4!1.0/(1.01 + 2.0*k*real(num_of_tdc)/real(L1))
  end subroutine getp
  
  subroutine get_lambda(num, lambda)
    implicit none
    integer, intent(in) :: num
    real(kind=8), intent(out) :: lambda
    
    if (num > L1) then
       lambda = lambda0*exp(-0.01*(num - L1))
    else
       lambda = lambda0
    end if
  end subroutine get_lambda
  
  subroutine output_to_file(pool, num, index)
    implicit none
    type (cell), intent(in) :: pool(:)
    integer, intent(in) :: index, num
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, num
       write(11, '(I2)', advance="no"), pool(i)%type
    end do
    write(11, '(2(F6.2))', advance="no"), p(0), t
    if (die_time < born_time) then 
       write(11, '(I8)', advance="no"), die_index
       write(11, '(I3)', advance="no"), -1
    else
       write(11, '(I8)', advance="no"), born_index
       write(11, '(I3)', advance="no"), 1
    end if
    close(11)
  end subroutine output_to_file

  subroutine cell_count(pool, num)
    implicit none
    type (cell), intent(in) :: pool(:)
    integer, intent(in) :: num
    integer i
    integer num_sc, num_tac, num_tdc

    num_sc = 0
    num_tac = 0
    num_tdc = 0
    do i = 1, num
       if ( pool(i)%type .eq. 1 ) then 
          num_sc = num_sc + 1
       elseif ( pool(i)%type .eq. 2 ) then 
          num_tac = num_tac + 1
       elseif ( pool(i)%type .eq. 3 ) then 
          num_tdc = num_tdc + 1
       else
          print *, 'woring'
          read(*,*)
       end if
    end do
    write(*, '(F10.2, 4I6, F10.4)'), t, num_sc, num_tac, num_tdc, num, p(1)
  end subroutine cell_count

  subroutine pool_print(unit, pool, num)
    implicit none
    integer, intent(in) :: unit
    type (cell), intent(in) :: pool(:)
    integer, intent(in) :: num
    
    integer i
    do i = 1, num
       write(unit, *), pool(i)%type
    end do
  end subroutine pool_print

  subroutine pool_init(pool, num)
    use random
    implicit none
    type (cell), intent(out) :: pool(:)
    integer, intent(out) :: num
    real(kind=8) u
    real(kind=8) r

    integer i
!!$    num = 20
!!$    pool(1:5)%type = 1
!!$    pool(6:10)%type = 2
!!$    pool(11:20)%type = 3
!!$    do i = 1, num
!!$       call ran2(u)
!!$       pool(i)%clock = u*lifespan
!!$    end do

    num = 10
    pool(1:10)%type = 1
    pool(1:10)%clock = lifespan(1)    
  end subroutine pool_init
  
end module paras
