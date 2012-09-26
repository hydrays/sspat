program FTCS1
  use setting
  implicit none
  open(unit=11,file='input.dat',status='unknown')
  open(unit=21,file='output.dat',status='unknown')
  open(unit=31,file='outputs/x.dat',status='unknown')
  open(unit=33,file='outputs/time.dat',status='unknown')
  open(unit=34,file='outputs/time2.dat',status='unknown')
  open(unit=35,file='outputs/t_trans.dat',status='unknown')

  CALL READ_IN          ! read-in data
  close (unit=11)

  CALL GRID             ! setting up the grid points in the domain

  CALL INITIAL_COND     ! setting up initial condition

  CALL BOUNDARY_COND    ! setting up boundary values

  CALL SOLVE            ! setting up system of equations

  CALL PRINTOUT         ! printing out the computational results

  CALL WRITE_OUT     ! write-out the basic data used for computation

  stop
end program FTCS1

