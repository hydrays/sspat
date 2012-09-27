program FTCS1
  use setting
  implicit none
  open(unit=11,file='input.dat',status='unknown')
  CALL READ_IN          ! read-in data
  close (unit=11)

  CALL GRID             ! setting up the grid points in the domain

  CALL INITIAL_COND     ! setting up initial condition

  CALL BOUNDARY_COND    ! setting up boundary values

  CALL SOLVE            ! setting up system of equations

  CALL WRITE_OUT     ! write-out the basic data used for computation

end program FTCS1

