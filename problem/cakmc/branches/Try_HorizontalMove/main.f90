program main
  use setting
  use kmcsample
  implicit none
  
  call read_xdata()

  if (useomp.eq.1) then
     print *, 'parallel not implemented yet.'
     read(*,*)
  else
     call kmc_sample_serial()
  end if

end program main
