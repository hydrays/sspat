module paras  
 
  integer, parameter :: NSpec = 5
  real(kind=8), parameter :: lifespan(NSpec) = (/1.0, 1.0, 1.0, 1.0, 1.0/)
  real(kind=8), parameter :: lambda0 = 10.0
  integer num, die_index, born_index
  integer num_of_dc, num_old
  real(kind=8) die_time, born_time
  real(kind=8) born_ruler
  real(kind=8) p(NSpec), v(NSpec)
  real(kind=8), parameter :: tend=1000
  real(kind=8), parameter :: L1=2000
  real(kind=8), parameter :: k1=0.9
  real(kind=8), parameter :: k2=8.0
  
  type cell
     integer type
     real(kind=8) clock
  end type cell

contains

  subroutine getp
    implicit none    
    p(1) = 1.0/(1.01 + k1*real(num_of_dc)/real(L1))
    p(2) = 0.55
    p(3) = 0.0
    p(4) = 0.6
    p(5) = 0.0
  end subroutine getp

  subroutine getv
    implicit none    
    v(1) = 2.0/(1.0 + k2*real(num_of_dc)/real(L1))
    v(2) = 1.0
    v(3) = 1.0
    v(4) = 1.0
    v(5) = 1.0
  end subroutine getv
  
  subroutine get_lambda(num, lambda)
    implicit none
    integer, intent(in) :: num
    real(kind=8), intent(out) :: lambda
    
    if (num > L1) then
       lambda = lambda0*exp(-0.01*(num - L1))
       !lambda = lambda0/(1.0+0.1*(num - L1))
    else
       lambda = lambda0
    end if
  end subroutine get_lambda
  
  subroutine output_to_file(pool, num, index, t)
    implicit none
    type (cell), intent(in) :: pool(:)
    integer, intent(in) :: index, num
    real(kind=8), intent(in) :: t
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

  subroutine cell_count(pool, num, t)
    implicit none
    type (cell), intent(in) :: pool(:)
    integer, intent(in) :: num
    real(kind=8), intent(in) :: t
    integer i, itype
    integer counter(NSpec)

    counter = 0
    do i = 1, num
       itype = pool(i)%type
       counter(itype) = counter(itype) + 1
    end do
    write(*, '(F10.2,6I6,2F10.4)'), t, counter, num, p(1),v(1)
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

    num = 2000
    pool(1:100)%type = 2
    pool(1:100)%clock = lifespan(1)
    pool(101:500)%type = 2
    pool(101:500)%clock = lifespan(2)
    pool(501:1990)%type = 3
    pool(501:1990)%clock = lifespan(3)
    pool(1991:2000)%type = 3
    pool(1991:2000)%clock = lifespan(4)
  end subroutine pool_init
  
end module paras
