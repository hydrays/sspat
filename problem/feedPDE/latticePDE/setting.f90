module setting
  implicit none
  integer, parameter :: nx = 256
  integer, parameter :: n = nx+2
  integer, parameter :: n1 = nx+3
  real x(0:n1)
  real phi_SC(0:n1), phi_PC(0:n1), phi_TC(0:n1), phi_MC(0:n1)
  real C1(0:n1), C2(0:n1), C3(0:n1), C4(0:n1)
  real p0(0:n1), pm(0:n1)
  real press(0:n1)
  real TGF(0:n1)
  real dx, dt, time
  real tw, te, qw1, qw2, qw3, qe1, qe2, qe3
  real bd10

  integer maxntimestp, ntimestp

  real xlen, Fo, tol,&
       alpha1, alpha2, alpha3, &
       alpha4, gain1, tfinal, tpinc, tm, &
       v_tc, gainm, l_d, p1, brange

  namelist /xdata/ xlen, Fo, tol, &
       alpha1, alpha2, alpha3, &
       alpha4, gain1, tfinal, tpinc, tm, &
       v_tc, gainm, l_d, p1, brange

contains

  subroutine update_TGF()
    implicit none
    integer i, k, shift_i
    
    do i = 2, n-1
       TGF(i) = 0.0
       do k = -brange, brange
          shift_i = k + i
          if ( shift_i .le. 1 ) then
             shift_i = shift_i + n - 2
          else if ( shift_i .ge. n ) then
             shift_i = shift_i - n + 2
          end if
          TGF(i) = TGF(i) + &
               bd10*phi_TC(shift_i)*exp(-real(abs(k))/brange)
       end do
       TGF(i) = TGF(i) / 100.0
    end do
    TGF(1) = TGF(n-1)
    TGF(n) = TGF(2)
  end subroutine update_TGF

  subroutine SOLVE
    implicit none
    real phi_SC_old(0:n1), phi_PC_old(0:n1)
    real phi_TC_old(0:n1), phi_MC_old(0:n1)
    integer i, j, file_index
    real dTdx_w, dTdx_e, dTdt, dTdx_west, dTdx_east
    real tp

    ! storing the values of temperature at nth time level
    phi_SC_old = phi_SC
    phi_PC_old = phi_PC
    phi_TC_old = phi_TC
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

       if (time .ge. tm) then
          phi_MC_old(nx/2) = 0.01
          tm = tm + tfinal
          !write(*,*), 'mutation time'
          !read(*,*)
       endif

       ! inclusion of source term
       call update_TGF()

       do i = 1, n
          press(i) = phi_SC_old(i) + phi_PC_old(i) + &
               phi_TC_old(i) + phi_MC_old(i)
       end do

       do i = 2, n-1
          p0(i) = 0.2 + 0.6 / (l_d + gain1*TGF(i))
          C1(i) = (2.0*p0(i)-1.0)*phi_SC_old(i)
          C2(i) = (2.0*(1-p0(i)))*phi_SC_old(i) + &
               (2.0*p1-1.0)*phi_PC_old(i)
          C3(i) = (2.0*(1-p1))*phi_PC_old(i) - &
               v_tc*phi_TC_old(i)
          C4(i) = 0.0
       enddo

       do i = 2, n-1
          dTdx_west = .5*(phi_SC_old(i) + phi_SC_old(i-1))*(press(i) - press(i-1))/dx
          dTdx_east = .5*(phi_SC_old(i+1) + phi_SC_old(i))*(press(i+1) - press(i))/dx  
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C1(i)
          phi_SC(i)    = phi_SC_old(i) + dTdt*dt

          dTdx_west = .5*(phi_PC_old(i) + phi_PC_old(i-1))*(press(i) - press(i-1))/dx
          dTdx_east = .5*(phi_PC_old(i+1) + phi_PC_old(i))*(press(i+1) - press(i))/dx  
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C2(i)
          phi_PC(i)    = phi_PC_old(i) + dTdt*dt

          dTdx_west = .5*(phi_TC_old(i) + phi_TC_old(i-1))*(press(i) - press(i-1))/dx
          dTdx_east = .5*(phi_TC_old(i+1) + phi_TC_old(i))*(press(i+1) - press(i))/dx  
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C3(i)
          phi_TC(i)    = phi_TC_old(i) + dTdt*dt

          dTdx_west = .5*(phi_MC_old(i) + phi_MC_old(i-1))*(press(i) - press(i-1))/dx
          dTdx_east = .5*(phi_MC_old(i+1) + phi_MC_old(i))*(press(i+1) - press(i))/dx
          dTdt      = alpha3*(dTdx_east - dTdx_west)/dx + C4(i)
          phi_MC(i)    = phi_MC_old(i) + dTdt*dt
       enddo
       ! boundary conditions

       call boundary_cond

       time = time + dt
       phi_SC_old = phi_SC
       phi_PC_old = phi_PC
       phi_TC_old = phi_TC
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

    do i=1,n
       write(11,'(10(e16.4e3))') phi_SC(i), phi_PC(i), &
            phi_TC(i), phi_MC(i), TGF(i)
    end do

    close(11)
  end subroutine output_to_file

  SUBROUTINE GRID
    implicit none
    integer i
    ! n-2 = number of control volumes the in x-direction
    dx   = xlen /(n-2)
    x(1) = 0
    do i = 2, n
       x(i) = x(i-1) + dx
    enddo

    open(10, file="out/x.dat")
    do i = 1, n
       write(10, '(f10.2)', advance="no"), x(i)
    enddo
    close(10)

  end SUBROUTINE GRID

  SUBROUTINE BOUNDARY_COND
    implicit none
    phi_SC(1) = phi_SC(n-1)
    phi_PC(1) = phi_PC(n-1)
    phi_TC(1) = phi_TC(n-1)
    phi_MC(1) = phi_MC(n-1)

    phi_SC(n) = phi_SC(2)
    phi_PC(n) = phi_PC(2)
    phi_TC(n) = phi_TC(2)
    phi_MC(n) = phi_MC(2)
  end SUBROUTINE BOUNDARY_COND

  SUBROUTINE INITIAL_COND
    implicit none
    integer i
    do i = 1, n
       phi_SC(i) = 0.1 + 0.05*sin(real(i))
       phi_PC(i) = 0.1
       phi_TC(i) = 0.2
       phi_MC(i) = 0.0
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

    bd10 = 10.0/real(brange)

    write(*, *), 'Control parameters...'
    write(*, '(a20, f10.2)'), 'xlen = ', xlen
    write(*, '(a20, f10.2)'), 'Fo = ', Fo
    write(*, '(a20, f10.2)'), 'tol = ', tol
    write(*, '(a20, f10.2)'), 'alpha1 = ', alpha1
    write(*, '(a20, f10.2)'), 'alpha2 = ', alpha2
    write(*, '(a20, f10.2)'), 'alpha3 = ', alpha3
    write(*, '(a20, f10.2)'), 'alpha4 = ', alpha4
    write(*, '(a20, f10.2)'), 'gain1 = ', gain1
    write(*, '(a20, f10.2)'), 'gainm = ', gainm
    write(*, '(a20, f10.2)'), 'tfinal = ', tfinal
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'p1 = ', p1
    write(*, '(a20, f10.2)'), 'tm = ', tm
    write(*, '(a20, f10.2)'), 'brange = ', brange

    open(9, file="out/parameter.dat")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    !write(9, '(a20, f10.4)'), 'xlen,', xlen
    !write(9, '(a20, f10.4)'), 'Fo,', Fo
    !write(9, '(a20, E10.4)'), 'tol,', tol
    write(9, '(a20, f10.4)'), 'alpha1,', alpha1
    write(9, '(a20, f10.4)'), 'alpha2,', alpha2
    write(9, '(a20, f10.4)'), 'alpha3,', alpha3
    write(9, '(a20, f10.4)'), 'alpha4,', alpha4
    write(9, '(a20, f10.4)'), 'gain1,', gain1
    write(9, '(a20, f10.4)'), 'gainm,', gainm
    write(9, '(a20, f10.4)'), 'tfinal,', tfinal
    write(9, '(a20, f10.4)'), 'tpinc,', tpinc
    write(9, '(a20, f10.4)'), 'tm,', tm
    write(9, '(a20, f10.4)'), 'v_tc,', v_tc
    write(9, '(a20, f10.4)'), 'p1,', p1
    write(9, '(a20, f10.4)'), 'l_d,', l_d
    write(9, '(a20, f10.4)'), 'brange,', brange

    close(8)
    close(9)
  end subroutine read_xdata

end module setting


