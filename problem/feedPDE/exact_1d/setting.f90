module setting
  implicit none
  integer boundary_type, maxntimestp, ntimestp, print_freq
  real, parameter :: pi = 3.14159265
  real alx, dx, tfinal, dt, time, alpha1, alpha2, Fo, tolstdy, tpinc
  integer, parameter :: nx = 1024
  integer, parameter :: n = nx+2
  integer, parameter :: n1 = nx+3
  real tw, te, qw1, qw2, qe1, qe2
  real x(0:n1)
  real phi1(0:n1)
  real phi2(0:n1)
  real C1(0:n1)
  real C2(0:n1)
  real p0(0:n1)
  real press(0:n1)
  real d(0:n1)
  real TGF(0:n1)
  real, parameter :: gain1 = 0.5
  real, parameter :: xi = 0.5

contains

  subroutine update_TGF(TGF(0:n1))
    implicit none
    real, intent(inout) :: TGF(0:n1)
    integer i
    real maxdif, dif
    real TGF_old(0:n1)

    TGF_old = TGF
    maxdif = huge(1.0)
    do while (maxdif > tolstdy)
       do i = 3, n-2
          dTdx_west = (TN1(i) - TN1(i-1))/dx
          dTdx_east = (TN1(i+1) - TN1(i))/dx
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C1(i)
          phi1(i)    = TN1(i) + dTdt*dt

          dTdx_west = (TN2(i) - TN2(i-1))/dx
          dTdx_east = (TN2(i+1) - TN2(i))/dx
          dTdt      = alpha2*(dTdx_east - dTdx_west)/dx + C2(i)
          phi2(i)    = TN2(i) + dTdt*dt
       enddo
       ! boundary conditions
       ! west
       i = 2
       if (boundary_type == 1) then
          print *, 'not implemented here!'
          stop
          dTdx_w = (-TN1(i+1) + 9*TN1(i) - 8*TN1(i-1))/(3*dx)
          dTdx_e = (TN1(i+1) - TN1(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = (-TN2(i+1) + 9*TN2(i) - 8*TN2(i-1))/(3*dx)
          dTdx_e = (TN2(i+1) - TN2(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       elseif (boundary_type == 2) then
          print *, 'not implemented here!'
          stop
       elseif (boundary_type == 3) then
          dTdx_w = qw1
          dTdx_e = (TN1(i+1) - TN1(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = qw2
          dTdx_e = (TN2(i+1) - TN2(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       end if
       ! east
       i = n-1
       if (boundary_type == 1) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = (8*TN(i+1) - 9*TN(i) + TN(i-1))/(3*dx)
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (boundary_type == 2) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = qe
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (boundary_type == 3) then
          dTdx_w = (TN1(i) - TN1(i-1))/dx
          dTdx_e = qe1
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = (TN2(i) - TN2(i-1))/dx
          dTdx_e = qe2
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       endif

       time = time + dt
       TN1 = phi1
       TN2 = phi2       
    end do
    
  end subroutine update_TGF
  
  subroutine SOLVE
    implicit none
    real TN1(0:n1), TN2(0:n1)
    integer i, j, file_index
    real dTdx_w, dTdx_e, dTdt, dTdx_west, dTdx_east
    real tp

    ! storing the values of temperature at nth time level
    TN1 = phi1
    TN2 = phi2

    time = 0.0
    tp = 0.0
    file_index = 0

    do ntimestp = 1, maxntimestp

       if (time .ge. tp) then
          call output_to_file(file_index)
          file_index = file_index + 1
          tp = tp + tpinc
          write(*, *), time
       endif

       ! inclusion of source term
       call update_TGF(TGF(0:n1))
       do i = 2, n-1
          p0(i) = 0.2 + 0.6 / (1.0 + gain1*TGF(i))
          press(i) = TN2(i) + TN1(i)
          d(i) = max(0.0, xi*(press(i) - 1))
          C1(i) = (2.0*p0(i)-1.0)*TN1(i) - d(i)*TN1(i)
          C2(i) = (2.0*(1-p0(i)))*TN1(i) - 0.2*TN2(i) - d(i)*TN2(i)
       enddo

       do i = 3, n-2
          dTdx_west = (TN1(i) - TN1(i-1))/dx
          dTdx_east = (TN1(i+1) - TN1(i))/dx
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C1(i)
          phi1(i)    = TN1(i) + dTdt*dt

          dTdx_west = (TN2(i) - TN2(i-1))/dx
          dTdx_east = (TN2(i+1) - TN2(i))/dx
          dTdt      = alpha2*(dTdx_east - dTdx_west)/dx + C2(i)
          phi2(i)    = TN2(i) + dTdt*dt
       enddo
       ! boundary conditions
       ! west
       i = 2
       if (boundary_type == 1) then
          print *, 'not implemented here!'
          stop
          dTdx_w = (-TN1(i+1) + 9*TN1(i) - 8*TN1(i-1))/(3*dx)
          dTdx_e = (TN1(i+1) - TN1(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = (-TN2(i+1) + 9*TN2(i) - 8*TN2(i-1))/(3*dx)
          dTdx_e = (TN2(i+1) - TN2(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       elseif (boundary_type == 2) then
          print *, 'not implemented here!'
          stop
       elseif (boundary_type == 3) then
          dTdx_w = qw1
          dTdx_e = (TN1(i+1) - TN1(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = qw2
          dTdx_e = (TN2(i+1) - TN2(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       end if
       ! east
       i = n-1
       if (boundary_type == 1) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = (8*TN(i+1) - 9*TN(i) + TN(i-1))/(3*dx)
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (boundary_type == 2) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = qe
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (boundary_type == 3) then
          dTdx_w = (TN1(i) - TN1(i-1))/dx
          dTdx_e = qe1
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi1(i) = TN1(i) + dTdt*dt

          dTdx_w = (TN2(i) - TN2(i-1))/dx
          dTdx_e = qe2
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi2(i) = TN2(i) + dTdt*dt
       endif

       time = time + dt
       TN1 = phi1
       TN2 = phi2

       !	calculation of heat fluxes at boundaries (flux is proportional to
       !   the negative of slope of the temperature curve)
       if (boundary_type == 1) then
          qw1 = -(-phi1(3) + 9*phi1(2) - 8*phi1(1)) /(3*dx)
          qe1 = -(8*phi1(n) - 9*phi1(n-1) + phi1(n-2)) /(3*dx)

          qw2 = -(-phi2(3) + 9*phi2(2) - 8*phi2(1)) /(3*dx)
          qe2 = -(8*phi2(n) - 9*phi2(n-1) + phi2(n-2)) /(3*dx)
       elseif (boundary_type == 2) then
          qw1 = -(-phi1(3) + 9*phi1(2) - 8*phi1(1)) /(3*dx)
          phi1(n) = (9*phi1(n-1) - phi1(n-2) + 3*qe1*dx) /8d0

          qw2 = -(-phi2(3) + 9*phi2(2) - 8*phi2(1)) /(3*dx)
          phi2(n) = (9*phi2(n-1) - phi2(n-2) + 3*qe2*dx) /8d0

       elseif (boundary_type == 3) then
          phi1(1) = (9*phi1(2) - phi1(3) + 3*qw1*dx) /8d0
          phi1(n) = (9*phi1(n-1) - phi1(n-2) + 3*qe1*dx) /8d0

          phi2(1) = (9*phi2(2) - phi2(3) + 3*qw2*dx) /8d0
          phi2(n) = (9*phi2(n-1) - phi2(n-2) + 3*qe2*dx) /8d0
       endif

    end do

  end subroutine SOLVE

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i=1,n
       write(11,'(3(e16.6))') phi1(i), phi2(i), p0(i)
    end do

    close(11)
  end subroutine output_to_file

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
    qw1 = 0.0 !alx/pi
    qw2 = 0.0 !alx/pi
    ! slope at west boundary (negative of the boundary flux)
    qe1 = 0.0 !-alx/pi
    qe2 = 0.0 !-alx/pi
    ! slope at east boundary (negative of the boundary flux)


    ! boundary_type: 1  - Dirichlet - all
    ! boundary_type: 2  - Dirichlet - left;   Neumann - right
    ! boundary_type: 3  - Neumann - right
    if (boundary_type == 1) then
       phi1(1) = tw
       phi1(n) = te
       phi2(1) = tw
       phi2(n) = te
    elseif (boundary_type == 2) then
       phi1(1) = tw
       phi2(1) = tw
    elseif (boundary_type == 3) then
       ! do nothing
    endif
  end SUBROUTINE BOUNDARY_COND

  SUBROUTINE INITIAL_COND
    implicit none
    integer i
    ! defining the initial conditions
    ! dsin(pi*x/alx)
    do i = 1, n
       phi1(i) = 0.1 + 0.05*sin(x(i))
       phi2(i) = 1.0
    enddo
  end SUBROUTINE INITIAL_COND

  subroutine READ_IN
    implicit none
    read(11,*) alx            ! actual length of domain
    read(11,*) boundary_type  ! flag for type of boundary conditions
    read(11,*) alpha1          ! thermal diffusivity
    read(11,*) alpha2          ! thermal diffusivity
    read(11,*) Fo             ! value of grid Fourier number for
    ! calculating the time step
    read(11,*) tolstdy        ! convergence criterion for steady-state
    read(11,*) tfinal         ! time at which the solution is desired
    read(11,*) tpinc     !
    ! calculation of time-step dt based on the stability condition
    dx = alx /(n-2)
    dt = Fo*dx**2 / max(alpha1, alpha2)

    maxntimestp = tfinal /dt
    dt = tfinal /maxntimestp   ! revised time-step
    Fo = max(alpha1, alpha2)*dt/dx**2        ! revised Fourier number
    if (tfinal < dt) then
       print*, 'Warning: final time is less than time-step!'
       stop
    endif
  end subroutine READ_IN

  subroutine WRITE_OUT
    implicit none

    open(unit=21,file='out/output.dat',status='unknown')

    if (boundary_type == 1) then
       write(21,*)'*Dirichlet bc is specified on all boundaries*'
    elseif (boundary_type == 2) then
       write(21,*)'*West-side Dirichlet, east-side Neumann conditions*'
    elseif (boundary_type == 3) then
       write(21,*)'*West-side Neumann, east-side Neumann conditions*'
    endif
    write(21,*) 'Thermal diffusivity', alpha1, alpha2
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
    write(21,*) 'The heat flux at West boundary is:', qw1, qw2
    write(21,*) 'The heat flux at East boundary is:', qe1, qe2
    close(21)
  end subroutine WRITE_OUT

end module setting


