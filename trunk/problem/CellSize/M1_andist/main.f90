program main
  use setting
  use cellgrow
  implicit none
  
  call read_xdata()
  
  call cell_grow()
  
end program main
