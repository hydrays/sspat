program main3d
  use setting
  use kmcsample3d
  implicit none

  call read_xdata()
  
  if (useomp.eq.1) then
     if (is64bit.eq.1) then
        call kmc_sample_omp3d64()
     else
        call kmc_sample_omp3d32()
        !call kmc_sample_omp32()
     end if
  else
     write(*, *), "Serial version of 3D not implemented yet...stop."
     stop
     !call kmc_sample_serial()
  end if

end program main3d
