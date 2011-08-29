program boosting
  use chem_data
  use random
  use nrtype
  implicit none

  real(kind=8) x(NSpec), x_flux(NSpec), dx(NSpec), mx(NReac), x_up(NSpec), x_down(NSpec)
  integer(I4B) i, j, k, is_nag
  real(kind=8) a(NReac), cuma_c(0:NReac), a_bar(NReac)
  real(kind=8) speed(NReac), kappa(NReac), highspeed
  integer(I4B) r(NReac), aorder(NReac), itemp
  logical cmask(NReac), rmask(NReac), rmask_old(NReac), RECORD_FLAG
  real(kind=8) tau, t, u, sum_a, dte, dt, dt2, mt, ct, load, loadtmp, mt2, pt
  real(kind=8), parameter :: EPS = 1e-20, epsilon = 0.1

  call ran_seed(sequence=1234)

  do k = 1, NSample
     t = 0.0
     kappa = 1.0
     x = xinit
     x_flux = 0
     dx = 0
     x_up = 0
     x_down = 0
     a_bar = 0
     pt = 0
     do j=1, NReac
        aorder(j) = j
     end do
!     mt = 100.0/highspeed
     mt = 0.2
     mt2 = mt
     ct = t + mt
     RECORD_FLAG = .FALSE.
     cmask = .true.
     rmask = .false.
     main_loop: do while (t < end_time)
        call ran1(u)
        sum_a = log(u)
        wait_cevent : do while (sum_a < -EPS)
           call getrate(x, a)
           a = kappa*a
           cmask = .true.
           mx = 1000000
           highspeed = EPS
           do j = 1, NReac
              do i = 1, NSpec
                 if (nu(i, j).ne.0) then
                    mx(j) = min(mx(j), x(i)/abs(nu(i, j)))
                 end if
              end do
              if (mx(j) > 10) then
                 cmask(j) = .false.
                 highspeed = max(highspeed, a(j)/(epsilon*mx(j)))
              end if
           end do
           dt2 = 1.0/highspeed
           r = 0
           cuma_c(0) = 0.0
           do j=1, NReac
              if (cmask(j)) then
                 cuma_c(j) = cuma_c(j-1) + a(j)
              else
                 cuma_c(j) = cuma_c(j-1)
              end if
           end do
           dte = -sum_a / cuma_c(NReac)
           dt = min(dte, dt2)
           sum_a = sum_a + cuma_c(NReac)*dt
           do j=1, NReac
              if (.not.cmask(j)) then
                 r(j) = poidev( a(j)*dt )
                 x = x + r(j)*nu(:,j)
                 if (rmask(j)) then
                    x_flux = x_flux + r(j)*abs(nu(:, j))
                    if (RECORD_FLAG) dx = dx + r(j)*nu(:, j)
                 end if
              end if
           end do
           if (RECORD_FLAG) then
              x_up = max(dx, x_up)
              x_down = min(dx, x_down)
           end if
           t = t + dt
           a_bar = a_bar + a*dt
           call checkx(x, is_nag)
           if (is_nag .eq. 1) then
              print *, 'nag1'
           end if
           if (t > pt*1.0) then
!              write(*, '(3(f20.5), 5(I8), 7(e12.4))'), t, mt, load, x, kappa
              pt = pt + 1
           end if
           if (t>end_time) then
              write(*, '(3(f20.5), 7(f8.1), 10(e12.4))'), t, mt, load, x, kappa
              exit main_loop
           end if
           if (t > ct) then
!! --- prepare cmask --- !!
              a_bar = a_bar/mt
              mx = 1000000
              do j = 1, NReac
                 do i = 1, NSpec
                    if (nu(i, j).ne.0) then
                       mx(j) = min(mx(j), x(i)/abs(nu(i, j)))
                    end if
                 end do
                 if (mx(j) > 10) then
                     a_bar(j) = a_bar(j)/(epsilon*mx(j))
                 end if
              end do
!! --- compute boost load --- !!
              if (RECORD_FLAG) then
                 dx = x_up - x_down
                 load = 0.0
                 do i = 1, NSpec
                    if (x_flux(i) .ne. 0) then
                       loadtmp = real(x_flux(i))/real(dx(i))
                       if(load.eq.0.0) then 
                          load = min(1000.0, loadtmp)
                       else if (loadtmp < load) then
                          load = loadtmp
                       end if
                    end if
                 end do
              end if
!! --- sort reactions --- !!
              do j=2, NReac
                 do i = j-1, 1, -1
                    if( a_bar(aorder(i)) < a_bar(aorder(i+1)) ) then
                       itemp = aorder(i)
                       aorder(i) = aorder(i+1)
                       aorder(i+1) = itemp
                    end if
                 end do
              end do
!! --- partition reactions --- !!
              rmask_old = rmask
              rmask=.false.
              rmask(aorder(1))=.true.
              loop_partition: do j=2, NReac
                 if (a_bar(aorder(j))>a_bar(aorder(j-1))*0.01) then
                    rmask(aorder(j)) = .true.
                 else
                    mt = 100.0/a_bar(aorder(j-1))
                    exit loop_partition
                 end if
              end do loop_partition
              if(rmask(aorder(NReac)) .or. a_bar(aorder(j))<0.001) then
!                 mt2 = 1.1*mt2
                 mt = mt2
!                 print *, 'no scale seperation'
!                 if (load < 10.0 .and. all(rmask.eqv.rmask_old)) then
                 if (load < 10.0) then
                    kappa = 1.0
                 end if
                 rmask = rmask_old
              else 
                 if (mt > 5.0*mt2) then
                    RECORD_FLAG = .FALSE.
                    mt2 = mt
!                    print *, 'check rmask', t, rmask
                 else
                    RECORD_FLAG = .true.
                    if (all(rmask.eqv.rmask_old)) then
                       mt2 = mt
                       if (load > 30.0 ) then
                          where(rmask)
                             kappa = 0.75*kappa
                          end where
                       else if (load < 10.0 ) then
                          kappa = 1.0
                       end if
                    end if
                 end if
              end if
!! --- reset all variables before exit --- !!
!              if (kappa(6).ne. kappa(5)) then
!              write(*, '(3(f20.5), 5(f8.1), 7(e12.4))'), t, mt, load, x, kappa
!              write(*, '(A, 73(f16.4))'), 'a_bar', a_bar
!              write(*, '(A, 73(f16.4))'), 'A_bar', a_bar(aorder)
!              write(*, '(A, 16(I12))'), 'x', x
!              write(*, '(A, 16(I12))'), 'flux', x_flux
!              write(*, '(A, 16(I12))'), 'dx', dx
!              write(*, '(A, 50(I6))'), 'order', aorder
!              print *, 'rmask', rmask
!              print *, 'cmask', cmask
!           end if
              ct = t + mt
              x_flux = 0
              dx = 0
              x_up = 0
              x_down = 0
              a_bar = 0
           end if
        end do wait_cevent
        call ran1(u)
        u = cuma_c(NReac)*u
        call assert( cuma_c(NReac).ne.0.0, 'cuma_c equal 0' )
        j = 1
        do while ( cuma_c(j) < u )
           j = j+1
        end do
        x = x + nu(:,j)
        if (rmask(j)) then
           x_flux = x_flux + abs(nu(:, j))
           if (RECORD_FLAG) then
              dx = dx + nu(:, j)
              x_up = max(dx, x_up)
              x_down = min(dx, x_down)
           end if
        end if
        call checkx(x, is_nag)
        if (is_nag .eq. 1) then
           print *, 'nag2',t 
        end if
     end do main_loop
  end do
end program boosting

subroutine checkx(x, is_nag)
  use chem_data
  use nrtype
  implicit none
  real(kind=8), intent(inout) :: x(NSpec)
  integer(I4B),  intent(out) :: is_nag
  is_nag = 0
  if(any(x <0) ) then
     is_nag = 1
     print *, 'nagtive!'
  end if
end subroutine checkx
