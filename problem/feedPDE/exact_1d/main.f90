program FTCS1
  use setting
  implicit none
  call read_xdata()

  CALL GRID             ! setting up the grid points in the domain

  CALL INITIAL_COND     ! setting up initial condition

  CALL BOUNDARY_COND    ! setting up boundary values

  CALL SOLVE            ! setting up system of equations

end program FTCS1

