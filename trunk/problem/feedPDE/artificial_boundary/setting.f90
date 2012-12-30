module setting
  implicit none
  integer, parameter :: nx = 512
  integer, parameter :: n = nx+2
  integer, parameter :: n1 = nx+3
  real x(0:n1)
  real phi_SC(0:n1), phi_TC(0:n1), phi_MC(0:n1)
  real C1(0:n1), C2(0:n1), C3(0:n1)
  real p0(0:n1), v0(0:n1)
  real press(0:n1)
  real q(0:n1)
  real d(0:n1)
  real TGF(0:n1)
  real dx, dt, time
  real tw, te, qw1, qw2, qw3, qe1, qe2, qe3

  integer btype, maxntimestp, ntimestp

  real xlen, Fo, tol,&
       alpha1, alpha2, alpha3, &
       alpha_TGF, gain1, xi, beta_prod, &
       beta_decay, tfinal, tpinc, v0max, v0min, tm, &
       v_tc, v_m, p_m, l_d

  namelist /xdata/ xlen, btype, Fo, tol, &
       alpha1, alpha2, alpha3, &
       alpha_TGF, gain1, xi, beta_prod, &
       beta_decay, tfinal, tpinc, v0max, v0min, tm, &
       v_tc, v_m, p_m, l_d

contains
  subroutine SOLVE
    implicit none
    real phi_SC_old(0:n1), phi_TC_old(0:n1), phi_MC_old(0:n1)
    integer i, j, file_index
    real dTdx_w, dTdx_e, dTdt, dTdx_west, dTdx_east
    real tp

    ! storing the values of temperature at nth time level
    phi_SC_old = phi_SC
    phi_MC_old = phi_MC

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

       do i = 2, n-1
          press(i) = phi_SC_old(i) + phi_MC_old(i)
!!$          if (press(i) > 0.65) then
!!$             press(1:i) = press(1:i) + 0.05
!!$          end if
          !press(i) = press(i) / (1.0 - 0.1*press(i))
          press(i) = press(i)
       end do

       do i = 2, n-1
          if (phi_MC_old(i) < 0.000000001) then
             phi_MC_old(i) = 0.0
          end if
       end do

       p0 = 0.55
       p0(150:200) = 0.5
       !q = 1.0
       q = 1.0 / (1.0 + exp(100.0*(press-1.0)))
       do i = 2, n-1
          if ( phi_MC_old(i) > 0.00001 .and. phi_MC_old(i+30) > 0.2) then
             p0(i) = 0.99
          end if
          !q(i) = 1.0
          !q(i) = 1.0/(1.0 + 2.0*press(i))
          d(i) = max(0.1, xi*(press(i)-1.0))
          v0(i) = 1.0
       end do

       do i = 2, n-1
          C1(i) = q(i)*v0(i)*(2.0*p0(i)-1.0)*phi_SC_old(i) - d(i)*phi_SC_old(i)
          C3(i) = q(i)*v_m*(2.0*p_m-1.0)*phi_MC_old(i) - d(i)*phi_MC_old(i)
       enddo

       do i = 3, n-2
          dTdx_west = (phi_SC_old(i) - phi_SC_old(i-1))/dx
          dTdx_east = (phi_SC_old(i+1) - phi_SC_old(i))/dx
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C1(i)
          phi_SC(i)    = phi_SC_old(i) + dTdt*dt

          dTdx_west = (phi_MC_old(i) - phi_MC_old(i-1))/dx
          dTdx_east = (phi_MC_old(i+1) - phi_MC_old(i))/dx
          dTdt      = alpha3*(dTdx_east - dTdx_west)/dx + C3(i)
          phi_MC(i)    = phi_MC_old(i) + dTdt*dt
       enddo
       ! boundary conditions
       ! west
       i = 2
       phi_SC(i) = 0.5
       phi_MC(i) = 0.0

       ! dTdx_w = qw1
       ! dTdx_e = (phi_SC_old(i+1) - phi_SC_old(i))/dx
       ! dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
       ! phi_SC(i) = phi_SC_old(i) + dTdt*dt
       
       ! dTdx_w = qw3
       ! dTdx_e = (phi_MC_old(i+1) - phi_MC_old(i))/dx
       ! dTdt   = alpha3*(dTdx_e - dTdx_w)/dx + C3(i)
       ! phi_MC(i) = phi_MC_old(i) + dTdt*dt
       ! east
       i = n-1
       phi_SC(i) = 0.0
       phi_MC(i) = 1.0
       ! dTdx_w = (phi_SC_old(i) - phi_SC_old(i-1))/dx
       ! dTdx_e = qe1
       ! dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
       ! phi_SC(i) = phi_SC_old(i) + dTdt*dt

       ! dTdx_w = (phi_MC_old(i) - phi_MC_old(i-1))/dx
       ! dTdx_e = qe3
       ! dTdt   = alpha3*(dTdx_e - dTdx_w)/dx + C3(i)
       ! phi_MC(i) = phi_MC_old(i) + dTdt*dt

       time = time + dt
       phi_SC_old = phi_SC
       phi_MC_old = phi_MC
    end do
  end subroutine SOLVE

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i=2,n-1
       !write(11,'(10(e16.4e3))') phi_SC(i), phi_TC(i), &
       !     phi_MC(i), p0(i), TGF(i), d(i), v0(i), press(i)
       write(11,'(10(e16.4e3))') phi_SC(i), &
            phi_MC(i), &
            q(i)*v0(i)*(2.0*p0(i)-1.0) - d(i), &
            p0(i), &
            q(i)*v_m - d(i), &
            press(i)
    end do

    close(11)
  end subroutine output_to_file

  SUBROUTINE GRID
    implicit none
    integer i
    ! n-2 = number of control volumes the in x-direction
    dx   = xlen /(n-2)
    x(1) = 0
    x(2) = x(1) + dx/2d0
    do i = 3, n-1
       x(i) = x(i-1) + dx
    enddo
    x(n) = x(n-1) + dx/2d0

    open(10, file="out/x.dat")
    do i = 2, n-1
       write(10, '(f10.2)', advance="no"), x(i)
    enddo
    close(10)

  end SUBROUTINE GRID

  SUBROUTINE BOUNDARY_COND
    implicit none
    ! Neumann-type
    qw1 = 0.0 !xlen/pi
    qw2 = 0.0 !xlen/pi
    qw3 = 0.0
    ! slope at west boundary (negative of the boundary flux)
    qe1 = 0.0 !-xlen/pi
    qe2 = 0.0 !-xlen/pi
    qe3 = 0.0
    ! slope at east boundary (negative of the boundary flux)

    !phi_SC(1) = 1.0
    !phi_SC(n) = 0.0

  end SUBROUTINE BOUNDARY_COND

  SUBROUTINE INITIAL_COND
    implicit none
    integer i, j 
    do i = 1, n/2
       phi_SC(i) = 0.5
       phi_MC(i) = 0.0
    enddo
    do i = n/2+1, n
       phi_SC(i) = 0.0
       phi_MC(i) = 1.0
    enddo
  end SUBROUTINE INITIAL_COND

  subroutine read_xdata()
    implicit none
    integer i

    open(8, file="input.dat", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    dx = xlen /(n-2)
    dt = Fo*dx**2 / max(alpha1, alpha2)

    maxntimestp = tfinal /dt
    dt = tfinal /maxntimestp   ! revised time-step
    Fo = max(alpha1, alpha2)*dt/dx**2        ! revised Fourier number
    if (tfinal < dt) then
       print*, 'Warning: final time is less than time-step!'
       stop
    endif

    write(*, *), 'Control parameters...'
    write(*, '(a20, f10.2)'), 'xlen = ', xlen
    write(*, '(a20, i10)'), 'btype = ', btype
    write(*, '(a20, f10.2)'), 'Fo = ', Fo
    write(*, '(a20, f10.2)'), 'tol = ', tol
    write(*, '(a20, f10.2)'), 'alpha1 = ', alpha1
    write(*, '(a20, f10.2)'), 'alpha2 = ', alpha2
    write(*, '(a20, f10.2)'), 'alpha3 = ', alpha3
    write(*, '(a20, f10.2)'), 'alpha_TGF = ', alpha_TGF
    write(*, '(a20, f10.2)'), 'gain1 = ', gain1
    write(*, '(a20, f10.2)'), 'xi = ', xi
    write(*, '(a20, f10.2)'), 'beta_prod = ', beta_prod
    write(*, '(a20, f10.2)'), 'beta_decay = ', beta_decay
    write(*, '(a20, f10.2)'), 'tfinal = ', tfinal
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'v0max = ', v0max
    write(*, '(a20, f10.2)'), 'tm = ', tm

    open(9, file="out/parameter.dat")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    !write(9, '(a20, f10.4)'), 'xlen,', xlen
    !write(9, '(a20, i10)'), 'btype,', btype
    !write(9, '(a20, f10.4)'), 'Fo,', Fo
    !write(9, '(a20, E10.4)'), 'tol,', tol
    write(9, '(a20, f10.4)'), 'alpha1,', alpha1
    write(9, '(a20, f10.4)'), 'alpha2,', alpha2
    write(9, '(a20, f10.4)'), 'alpha3,', alpha3
    write(9, '(a20, f10.4)'), 'alpha_TGF,', alpha_TGF
    write(9, '(a20, f10.4)'), 'gain1,', gain1
    write(9, '(a20, f10.4)'), 'xi,', xi
    write(9, '(a20, f10.4)'), 'beta_prod,', beta_prod
    write(9, '(a20, f10.4)'), 'beta_decay,', beta_decay
    write(9, '(a20, f10.4)'), 'tfinal,', tfinal
    write(9, '(a20, f10.4)'), 'tpinc,', tpinc
    write(9, '(a20, f10.4)'), 'v0max,', v0max
    write(9, '(a20, f10.4)'), 'v0min,', v0min
    write(9, '(a20, f10.4)'), 'tm,', tm
    write(9, '(a20, f10.4)'), 'v_tc,', v_tc
    write(9, '(a20, f10.4)'), 'v_m,', v_m
    write(9, '(a20, f10.4)'), 'p_m,', p_m
    write(9, '(a20, f10.4)'), 'l_d,', l_d

    close(8)
    close(9)
  end subroutine read_xdata

end module setting


