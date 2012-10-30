module setting
  implicit none
  integer, parameter :: nx = 256
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

  subroutine update_TGF()
    implicit none
    integer i
    real maxdif_TGF, dif_TGF, dt_TGF
    real TGF_old(0:n1)
    real C_TGF(0:n1)
    real dTdx_w, dTdx_e, dTdt, dTdx_west, dTdx_east
    integer i_counter
    integer boundary_condition_TGF

    TGF_old = TGF
    maxdif_TGF = huge(1.0)
    dt_TGF = Fo*dx**2 / alpha_TGF
    i_counter = 0

    do while (maxdif_TGF > 1e-5)

    do i = 2, n-1
       C_TGF(i) = beta_prod*phi_TC(i) - beta_decay*TGF(i)
    end do

       do i = 3, n-2
          dTdx_west = (TGF_old(i) - TGF_old(i-1))/dx
          dTdx_east = (TGF_old(i+1) - TGF_old(i))/dx
          dTdt      = alpha_TGF*(dTdx_east - dTdx_west)/dx + C_TGF(i)
          TGF(i)   = TGF_old(i) + dTdt*dt_TGF
       enddo
       ! boundary conditions
       ! west Direchlet
       i = 2
       !dTdx_w = (-TGF_old(i+1) + 9*TGF_old(i) - 8.0*TGF_old(i-1))/(3*dx)
       dTdx_w = 0.0
       dTdx_e = (TGF_old(i+1) - TGF_old(i))/dx
       dTdt   = alpha_TGF*(dTdx_e - dTdx_w)/dx + C_TGF(i)
       TGF(i) = TGF_old(i) + dTdt*dt_TGF
       ! east Direchlet
       i = n-1
       dTdx_w = (TGF_old(i) - TGF(i-1))/dx
       !dTdx_e = (8*TGF_old(i+1) - 9*TGF_old(i) + TGF_old(i-1))/(3*dx)
       dTdx_e = 0.0
       dTdt   = alpha_TGF*(dTdx_e - dTdx_w)/dx + C_TGF(i)
       TGF(i) = TGF_old(i) + dTdt*dt_TGF

       maxdif_TGF = 0.0
       do i = 2, n-1
          dif_TGF = TGF(i) - TGF_old(i)
          if (abs(dif_TGF) > maxdif_TGF) then
             maxdif_TGF = abs(dif_TGF)
          endif
       enddo
       TGF_old = TGF
       i_counter = i_counter + 1
       !write(*,*), time, i_counter, maxdif_TGF
       !write(*,'(6(f10.6))'), C_TGF(100:105)
       !write(*,'(6(f10.6))'), TGF(100:105)
    end do
    !write(*,*), i_counter, maxdif_TGF
  end subroutine update_TGF

  subroutine SOLVE
    implicit none
    real phi_SC_old(0:n1), phi_TC_old(0:n1), phi_MC_old(0:n1)
    integer i, j, file_index
    real dTdx_w, dTdx_e, dTdt, dTdx_west, dTdx_east
    real tp

    ! storing the values of temperature at nth time level
    phi_SC_old = phi_SC
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

       do i = 2, n-1
          !press(i) = phi_TC_old(i) + phi_SC_old(i) + phi_MC_old(i)
          press(i) = phi_SC_old(i) + phi_MC_old(i)
!!$          if (press(i) > 0.65) then
!!$             press(1:i) = press(1:i) + 0.05
!!$          end if
          !press(i) = press(i) / (1.0 - 0.1*press(i))
          press(i) = press(i)
       end do

       do i = 2, n-1
          if (phi_MC_old(i) < tol) then
             phi_MC_old(i) = 0.0
          end if
       end do

       do i = 2, n-1
          p0(i) = 1.0 / (l_d + (gain1*TGF(i))**1.0)
!!$          if ( press(i) > 0.8) then
!!$             q(i) = 0.0
!!$          else
!!$             q(i) = 1.0
!!$          end if
          q(i) = 1.0 / (1.0 + exp(100.0*(press(i)-0.8)))
          !q(i) = 1.0/(1.0 + 2.0*press(i))
          d(i) = max(0.0, xi*(press(i)-0.8))
          v0(i) = 1 / (1.0/v0max + gain1*TGF(i)*(1.0/v0min - 1.0/v0max))
          C1(i) = q(i)*v0(i)*(2.0*p0(i)-1.0)*phi_SC_old(i) - d(i)*phi_SC_old(i)
          C2(i) = (2.0*(1-p0(i)))*phi_SC_old(i) - &
               v_tc*phi_TC_old(i)
!!$          if ( time > 15 ) then
!!$             if ( sum(phi_MC_old(2:n-1)) > 0.01 ) then
!!$                !print *, time, sum(phi_MC_old(2:n-1))
!!$                C2(2:2) = C2(2:2) - 10.0*phi_TC_old(2:2)
!!$             else
!!$                phi_MC_old(2:n-1) = 0.0
!!$             end if
!!$          end if
          C3(i) = q(i)*v_m*(2.0*p_m-1.0)*phi_MC_old(i) - d(i)*phi_MC_old(i)
       enddo

       do i = 2, n-1
          if (time > 120.0) then
             C2(i) = C2(i) - 1000.0*phi_MC_old(i)*phi_TC_old(i)
          end if
       end do

       do i = 3, n-2
          dTdx_west = (phi_SC_old(i) - phi_SC_old(i-1))/dx
          dTdx_east = (phi_SC_old(i+1) - phi_SC_old(i))/dx
          dTdt      = alpha1*(dTdx_east - dTdx_west)/dx + C1(i)
          phi_SC(i)    = phi_SC_old(i) + dTdt*dt

          dTdx_west = (phi_TC_old(i) - phi_TC_old(i-1))/dx
          dTdx_east = (phi_TC_old(i+1) - phi_TC_old(i))/dx
          dTdt      = alpha2*(dTdx_east - dTdx_west)/dx + C2(i)
          phi_TC(i)    = phi_TC_old(i) + dTdt*dt

          dTdx_west = (phi_MC_old(i) - phi_MC_old(i-1))/dx
          dTdx_east = (phi_MC_old(i+1) - phi_MC_old(i))/dx
          dTdt      = alpha3*(dTdx_east - dTdx_west)/dx + C3(i)
          phi_MC(i)    = phi_MC_old(i) + dTdt*dt
       enddo
       ! boundary conditions
       ! west
       i = 2
       if (btype == 1) then
          print *, 'not implemented here!'
          stop
          dTdx_w = (-phi_SC_old(i+1) + 9*phi_SC_old(i) -&
               8*phi_SC_old(i-1))/(3*dx)
          dTdx_e = (phi_SC_old(i+1) - phi_SC_old(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi_SC(i) = phi_SC_old(i) + dTdt*dt

          dTdx_w = (-phi_TC_old(i+1) + 9*phi_TC_old(i) -&
               8*phi_TC_old(i-1))/(3*dx)
          dTdx_e = (phi_TC_old(i+1) - phi_TC_old(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi_TC(i) = phi_TC_old(i) + dTdt*dt
       elseif (btype == 2) then
          print *, 'not implemented here!'
          stop
       elseif (btype == 3) then
          dTdx_w = qw1
          dTdx_e = (phi_SC_old(i+1) - phi_SC_old(i))/dx
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi_SC(i) = phi_SC_old(i) + dTdt*dt

          dTdx_w = qw2
          dTdx_e = (phi_TC_old(i+1) - phi_TC_old(i))/dx
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi_TC(i) = phi_TC_old(i) + dTdt*dt

          dTdx_w = qw3
          dTdx_e = (phi_MC_old(i+1) - phi_MC_old(i))/dx
          dTdt   = alpha3*(dTdx_e - dTdx_w)/dx + C3(i)
          phi_MC(i) = phi_MC_old(i) + dTdt*dt
       end if
       ! east
       i = n-1
       if (btype == 1) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = (8*TN(i+1) - 9*TN(i) + TN(i-1))/(3*dx)
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (btype == 2) then
          print *, 'not implemented here!'
          stop
          !dTdx_w = (TN(i) - TN(i-1))/dx
          !dTdx_e = qe
          !dTdt   = alpha*(dTdx_e - dTdx_w)/dx + C(i)
          !phi(i) = TN(i) + dTdt*dt
       elseif (btype == 3) then
          dTdx_w = (phi_SC_old(i) - phi_SC_old(i-1))/dx
          dTdx_e = qe1
          dTdt   = alpha1*(dTdx_e - dTdx_w)/dx + C1(i)
          phi_SC(i) = phi_SC_old(i) + dTdt*dt

          dTdx_w = (phi_TC_old(i) - phi_TC_old(i-1))/dx
          dTdx_e = qe2
          dTdt   = alpha2*(dTdx_e - dTdx_w)/dx + C2(i)
          phi_TC(i) = phi_TC_old(i) + dTdt*dt

          dTdx_w = (phi_MC_old(i) - phi_MC_old(i-1))/dx
          dTdx_e = qe3
          dTdt   = alpha3*(dTdx_e - dTdx_w)/dx + C3(i)
          phi_MC(i) = phi_MC_old(i) + dTdt*dt
       endif

       time = time + dt
       phi_SC_old = phi_SC
       phi_TC_old = phi_TC
       phi_MC_old = phi_MC

       ! calculation of heat fluxes at boundaries (flux is proportional to
       ! the negative of slope of the temperature curve)
       if (btype == 1) then
          qw1 = -(-phi_SC(3) + 9*phi_SC(2) - 8*phi_SC(1)) /(3*dx)
          qe1 = -(8*phi_SC(n) - 9*phi_SC(n-1) + phi_SC(n-2)) /(3*dx)

          qw2 = -(-phi_TC(3) + 9*phi_TC(2) - 8*phi_TC(1)) /(3*dx)
          qe2 = -(8*phi_TC(n) - 9*phi_TC(n-1) + phi_TC(n-2)) /(3*dx)
       elseif (btype == 2) then
          qw1 = -(-phi_SC(3) + 9*phi_SC(2) - 8*phi_SC(1)) /(3*dx)
          phi_SC(n) = (9*phi_SC(n-1) - phi_SC(n-2) + 3*qe1*dx) /8d0

          qw2 = -(-phi_TC(3) + 9*phi_TC(2) - 8*phi_TC(1)) /(3*dx)
          phi_TC(n) = (9*phi_TC(n-1) - phi_TC(n-2) + 3*qe2*dx) /8d0

       elseif (btype == 3) then
          phi_SC(1) = (9*phi_SC(2) - phi_SC(3) + 3*qw1*dx) /8d0
          phi_SC(n) = (9*phi_SC(n-1) - phi_SC(n-2) + 3*qe1*dx) /8d0

          phi_TC(1) = (9*phi_TC(2) - phi_TC(3) + 3*qw2*dx) /8d0
          phi_TC(n) = (9*phi_TC(n-1) - phi_TC(n-2) + 3*qe2*dx) /8d0

          phi_MC(1) = (9*phi_MC(2) - phi_MC(3) + 3*qw3*dx) /8d0
          phi_MC(n) = (9*phi_MC(n-1) - phi_MC(n-2) + 3*qe3*dx) /8d0
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

    do i=2,n-1
       !write(11,'(10(e16.4e3))') phi_SC(i), phi_TC(i), &
       !     phi_MC(i), p0(i), TGF(i), d(i), v0(i), press(i)
       write(11,'(10(e16.4e3))') phi_SC(i), phi_TC(i), &
            phi_MC(i), q(i)*v0(i), q(i)*v_m, &
            !q(i)*v0(i)*(2.0*p0(i)-1.0), &
            p0(i), &
            q(i)*v_m*(2.0*p_m-1.0), press(i)
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

    ! defining the boundary conditions
    tw = 1

    ! temperature at west (left) boundary
    te = 0
    ! temperature at east (right) boundary

    ! Neumann-type
    qw1 = 0.0 !xlen/pi
    qw2 = 0.0 !xlen/pi
    qw3 = 0.0
    ! slope at west boundary (negative of the boundary flux)
    qe1 = 0.0 !-xlen/pi
    qe2 = 0.0 !-xlen/pi
    qe3 = 0.0
    ! slope at east boundary (negative of the boundary flux)


    ! btype: 1  - Dirichlet - all
    ! btype: 2  - Dirichlet - left;   Neumann - right
    ! btype: 3  - Neumann - right
    if (btype == 1) then
       phi_SC(1) = tw
       phi_SC(n) = te
       phi_TC(1) = tw
       phi_TC(n) = te
    elseif (btype == 2) then
       phi_SC(1) = tw
       phi_TC(1) = tw
    elseif (btype == 3) then
       ! do nothing
    endif

    TGF(n) = 0.0
    TGF(1) = 0.0
  end SUBROUTINE BOUNDARY_COND

  SUBROUTINE INITIAL_COND
    implicit none
    integer i
    ! defining the initial conditions
    ! dsin(pi*x/xlen)
    do i = 1, n
       phi_SC(i) = 0.1 + 0.05*sin(x(i))
       phi_TC(i) = 1.0
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


