module setting
  implicit none
  integer boundary_type, maxntimestp, ntimestp, print_freq
  real, parameter :: pi = 3.14159265
  real alx, dx, tfinal, dt, time, alpha, Fo, tolstdy
  integer, parameter :: nx = 64
  integer, parameter :: n = nx+2
  integer, parameter :: n1 = nx+3
  real tw, te, qw, qe 
  real x(0:n1)
  real phi(0:n1)
  real C(0:n1)

contains
subroutine READ_IN
  implicit none
  read(11,*) alx            ! actual length of domain
  read(11,*) boundary_type  ! flag for type of boundary conditions
  read(11,*) alpha          ! thermal diffusivity
  read(11,*) Fo             ! value of grid Fourier number for 
  ! calculating the time step
  read(11,*) tolstdy        ! convergence criterion for steady-state
  read(11,*) tfinal         ! time at which the solution is desired
  read(11,*) print_freq     ! 
  ! calculation of time-step dt based on the stability condition
  dx = alx /(n-2)
  dt = Fo*dx**2 /alpha
  maxntimestp = tfinal /dt
  dt = tfinal /maxntimestp   ! revised time-step
  Fo = alpha*dt/dx**2        ! revised Fourier number
  if (tfinal < dt) then
     print*, 'Warning: final time is less than time-step!'
     stop
  endif
end subroutine READ_IN

subroutine WRITE_OUT
  implicit none

  if (boundary_type == 1) then
     write(21,*)'*Dirichlet bc is specified on all boundaries*'
  elseif (boundary_type == 2) then
     write(21,*)'*West-side Dirichlet, east-side Neumann conditions*'
  endif
  write(21,*) 'Thermal diffusivity', alpha
  write(21,*) 'Grid Fourier number, Fo =', Fo
  write(21,*) 'Length of the domain', alx
  write(21,*) 'Number of control volumes:', nx
  write(21,*) 'dx =', dx
  write(21,*) 'dt =', dt
  write(21,*) 'No. of Time step = ', ntimestp
  write(21,*) 'Final time = ', time
  write(21,*) 'convergence criterion for steady-state solution:', tolstdy
  write(21,*) ''
  write(21,*) 'The temperature at West boundary is:', tw
  write(21,*) 'The temperature at East boundary is:', te
  write(21,*) 'The heat flux at West boundary is:', qw
  write(21,*) 'The heat flux at East boundary is:', qe
  return
end subroutine WRITE_OUT

SUBROUTINE GRID
  implicit none
  integer i
  ! n-2 = number of control volumes the in x-direction
  dx   = alx /(n-2)
  x(1) = 0
  x(2) = x(1) + dx/2d0
  do i = 3, n-1
     x(i) = x(i-1) + dx
  enddo
  x(n) = x(n-1) + dx/2d0
  return
end SUBROUTINE GRID

SUBROUTINE BOUNDARY_COND
  implicit none
  
  ! defining the boundary conditions
  tw = 1		
  ! temperature at west (left) boundary 
  te = 0		
  ! temperature at east (right) boundary
  
  ! Neumann-type
  qw = alx/pi	
  ! slope at west boundary (negative of the boundary flux)
  qe = 0 !-alx/pi	
  ! slope at east boundary (negative of the boundary flux)


  ! boundary_type: 1  - Dirichlet - all
  ! boundary_type: 2  - Dirichlet - left;   Neumann - right
  if (boundary_type == 1) then
     phi(1) = tw
     phi(n) = te
  elseif (boundary_type == 2) then
     phi(1) = tw
  endif
end SUBROUTINE BOUNDARY_COND

SUBROUTINE INITIAL_COND
  implicit none
  integer i
  ! defining the initial conditions
  ! dsin(pi*x/alx)
  do i = 2, n-1
     phi(i) = 0.0
  enddo
end SUBROUTINE INITIAL_COND

subroutine SOLVE
  implicit none
  real TN(0:n1), maxdif
  integer monitor, jcounter, i, j
  real dTdx_w, dTdx_west, dTdx_e, dTdx_east, dif, dTdt

  ! storing the values of temperature at nth time level
  TN = phi
  ! inclusion of source term
  do i = 2, n-1
     C(i) = 0.0
  enddo
  
  ntimestp = 0
  maxdif = 1d9
  time = 0d0
  write(33,101) time        ! continuous time
  write(34,101) time        ! periodic time
  monitor = nx/4
  write(35,*) phi(monitor)
  
10 if ((maxdif > tolstdy).and.(ntimestp < maxntimestp)) then
     do i = 3, n-2
        dTdx_west = (TN(i) - TN(i-1))/dx
        dTdx_east = (TN(i+1) - TN(i))/dx
        dTdt      = alpha*(dTdx_east - dTdx_west)/dx + C(i)
        phi(i)    = TN(i) + dTdt*dt
     enddo
     !	   boundary conditions
     !	   west
     i = 2
     dTdx_w = (-TN(i+1) + 9*TN(i) - 8*TN(i-1))/(3*dx)
     dTdx_e = (TN(i+1) - TN(i))/dx
     dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
     phi(i) = TN(i) + dTdt*dt
     !	   east
     i = n-1
     if (boundary_type == 1) then
        dTdx_w = (TN(i) - TN(i-1))/dx
        dTdx_e = (8*TN(i+1) - 9*TN(i) + TN(i-1))/(3*dx)
        dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
        phi(i) = TN(i) + dTdt*dt
     elseif (boundary_type == 2) then
        dTdx_w = (TN(i) - TN(i-1))/dx
        dTdx_e = qe
        dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
        phi(i) = TN(i) + dTdt*dt
     endif

     time     = time + dt
     ntimestp = ntimestp + 1
     
     maxdif = 0
     do i = 2, n-1
        dif = phi(i) - TN(i)
        if (abs(dif) > maxdif) then
           maxdif = abs(dif)
        endif
     enddo

     !   updating the temperature
     TN = phi
     
     write(33,101) time              ! continuous time
     write(6,102) ntimestp, maxdif

     !   printing out transient temperature at selected point of domain
     j = jcounter + 1
     jcounter = mod(j,print_freq)
     if (jcounter == 0) then
        write(34,101) (time + dt)    ! periodic time
        write(35,*) phi(monitor)
     endif
     go to 10
  endif

  !	calculation of heat fluxes at boundaries (flux is proportional to 
  !                      the negative of slope of the temperature curve)
  if (boundary_type == 1) then
     qw = -(-phi(3) + 9*phi(2) - 8*phi(1)) /(3*dx)
     qe = -(8*phi(n) - 9*phi(n-1) + phi(n-2)) /(3*dx)
  elseif (boundary_type == 2) then
     qw = -(-phi(3) + 9*phi(2) - 8*phi(1)) /(3*dx)
     phi(n) = (9*phi(n-1) - phi(n-2) + 3*qe*dx) /8d0
  endif
  
101 format(e14.7)
102 format(i7,5x,e14.7)
  
  return
end subroutine SOLVE

SUBROUTINE PRINTOUT
  implicit none
  integer i
  do i=1,n
     write(31,'(2(e13.6))') x(i), phi(i)
  end do
end SUBROUTINE PRINTOUT

end module setting


