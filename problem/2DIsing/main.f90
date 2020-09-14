program main
  use setting
  implicit none

  real t, tp, u
  integer output_index
  integer i, j
  integer, allocatable :: seed(:)
  real E0, E1
  integer size
  
  call read_xdata()

  call random_seed(size=size)
  allocate(seed(size))
  seed(:) = iseed
  call random_seed(put=seed)

  open (unit = 100, file='./out/logfile', action="write")

  call init()
  call padding()
  
  t = 0.0
  tp = 0.0
  output_index = 0

  do while (t < tend)
     if (t .ge. tp) then
        !call cell_stat(t)
        call output_to_file(output_index)
        output_index = output_index + 1
        tp = tp + tpinc
     end if

     call get_energy(E0)
     
     ! source: j = n + FLOOR((m+1-n)*u)
     call random_number(u)
     i = pad+1 + FLOOR((Lbox-2*pad)*u)
     call random_number(u)
     j = pad+1 + FLOOR((Lbox-2*pad)*u)

     cmat(i, j) = -1 * cmat(i, j)
     call get_energy(E1)

     call random_number(u)     
     if ( u < prob(E0, E1) ) then
        ! do nothing
        call padding()
     else
        cmat(i, j) = -1 * cmat(i, j)
     end if

     t = t + dt
  end do

  close(unit=100)

end program main
