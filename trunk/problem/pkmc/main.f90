program main
  use omp_lib
  use nrtype
  use random
  use setting
  implicit none

  real t, tau, tp, tm, u
  integer output_index, i, j, active_index, temp_num
  real tau_temp
  integer k
  open (unit = 100, file='./out/logfile', action="write")
  call ran_seed(sequence=12345)

  cmat(1:L, 1)%type = 1
  cmat(1:L, 2:3)%type = 2
  cmat(1:L, 3:4)%type = 3
!  cmat(1:L, 5:50)%type = 1
  do i = 1, L
     call ran2(u)
     cmat(i, 1)%gene1 = 0.5*u
  end do

  t = 0.0
  tp = 0.0
  tm = 100000.0
  output_index = 0
  NP = 0.0
  NT = 0.0
  do i = 1, L
     call expdev(u)
     NP(i) = u
  end do
  cmat(0, :) = cmat(L, :)
  cmat(L+1, :) = cmat(1, :)

  npack = 0
  do i = 0, L+1
     do j = 1, H
        if ( cmat(i,j)%type .ne. 0 ) then
           npack(i) = npack(i) + 1
        end if
     end do
  end do
     
  TGFbeta = 0.
  do i = 1, L
     temp_num = 0.
     do j = 1, H
        if ( cmat(i,j)%type .eq. 3 ) then
           temp_num = temp_num + 1.
        end if
     end do
     !temp_num = 10.0*temp_num / b
     if (temp_num > 0) then
        do j = 0, b
           if ( j .eq. 0) then
              TGFbeta(i) = TGFbeta(i) + temp_num
           else
              TGFbeta(i+j) = TGFbeta(i+j) + temp_num*exp(-real(j)/b)
              TGFbeta(i-j) = TGFbeta(i-j) + temp_num*exp(-real(j)/b)
           end if
        end do
     end if
  end do
  do j = 1, b
     TGFbeta(j) = TGFbeta(j) + TGFbeta(L+j)
     TGFbeta(L-j+1) = TGFbeta(L-j+1) + TGFbeta(1-j)
  end do
  TGFbeta(-b:0) = TGFbeta(L-b:L)
  TGFbeta(L+1:L+b+1) = TGFbeta(1:b+1)

  D_TGFbeta = 0.0
  do j = 1, 1+2*b
     D_TGFbeta(j) = exp(-real(abs(j-b-1))/b)
  end do
  
  do while (t < tend)
     if (t .ge. tp) then
        call output_to_file(output_index)
        call cell_stat(t)
!        print *, t
        output_index = output_index + 1
        tp = tp + 5.0
     end if

     if (t .ge. tm) then
        cmat(100:105, 1:10)%type = 4
        tm = tm + 600000.0
     end if

     tau = 11111111111.1
     k = 0
     do i = 1, L
        vr = exp(real(npack(i)-npack(i+1)))
        vl = exp(real(npack(i)-npack(i-1)))
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
