program main
  use omp_lib
  use nrtype
  use random
  use setting
  implicit none

  real t, tau, tp, tm, u
  integer output_index, i, j, active_index
  real tau_temp
  integer k
  open (unit = 100, file='./out/logfile', action="write")
  call ran_seed(sequence=12345)

  call init_cell_pool()

  t = 0.0
  tp = 0.0
  tm = 100000.0
  output_index = 0
  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(output_index)
        call cell_stat(t)
!        print *, t
        output_index = output_index + 1
        tp = tp + 1.0
     end if

     if (t .ge. tm) then
        cmat(100:105, 1:10)%type = 4
        tm = tm + 600000.0
     end if

     tau = 11111111111.1
     k = 0
     do i = 1, L
        !vr = exp(real(npack(i)-npack(i+1)))
        !vl = exp(real(npack(i)-npack(i-1)))
        vr = max(0.0, 100.0*real(npack(i)-npack(i+1)))
        vl = max(0.0, 100.0*real(npack(i)-npack(i-1)))
        a(i) = vr + vl + npack(i)*v
        tau_temp = ( NP(i) - NT(i) ) / a(i)
        if ( tau_temp < tau) then
           tau = tau_temp
           k = i
        end if
     end do
     if ( k .le. 0 ) then
        write(*,*), 'error'
        read(*,*)
     end if
     do i = 1, L
        NT(i) = NT(i) + a(i)*tau
     end do
     call expdev(u)
     NP(k) = NP(k) + u
     t = t + tau
     do i = 1, L
        if ( i .eq. k ) then
           call cell_event(k)
        elseif ( abs(i - k) > 100 ) then
           call cell_restack(i)
        end if
     end do
     ! perodic boundary condition
     if ( k .eq. 1 ) then
        cmat(L, :) = cmat(0, :)
        cmat(L+1, :) = cmat(1, :)
        npack(L) = npack(0)
        npack(L+1) = npack(1)
     end if
     if ( k .eq. 2 ) then
        cmat(L+1, :) = cmat(1, :)
        npack(L+1) = npack(1)
     end if
     if ( k .eq. L ) then
        cmat(1, :) = cmat(L+1, :)
        cmat(0, :) = cmat(L, :)
        npack(1) = npack(L+1)
        npack(0) = npack(L)
     end if
     if ( k .eq. L-1 ) then
        cmat(0, :) = cmat(L, :)
        npack(0) = npack(L)
     end if
     if ( k .le. 2*b+1+1 ) then
        TGFbeta(L-b:L+b+1) = TGFbeta(-b:b+1)
     end if
     if ( k .ge. L - 2*b-1 ) then
        TGFbeta(-b:b+1) = TGFbeta(L-b:L+b+1)
     end if
  end do

  close(unit=100)
end program main
