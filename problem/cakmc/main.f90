program main
  use setting
  use kmcsample
  implicit none
  
  call read_xdata()

  if (useomp.eq.1) then
     if (is64bit.eq.1) then
        call kmc_sample_omp64()
     else
        call kmc_sample_omp32()
     end if
  else
     call kmc_sample_serial()
  end if

end program main
