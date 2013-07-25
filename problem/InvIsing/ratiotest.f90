program Ratiotest
  implicit none
  integer, parameter :: p = 4 !/*system size*/
  integer, parameter :: n = 20000	!/*number of samples*/
  real, parameter :: B = 1.0 !/*inverse temperature*/
  integer i,j,k,l
  real J1(p,p),h1(p), J2(p,p),h2(p), J3(p,p),h3(p)
  real samples(n, p)
  integer x(p)
  real Z1, Z2, Z3, f
  real R, p1, p2

  open(unit = 21, file="sample.dat", action = "read")
  do l = 1, n
     read(21, *) samples(l, :)
  end do
  close(21)

  open(file="parameters.dat",unit=41,action="read")
  open(file="parametersest.dat",unit=42,action="read")
  do j = 1, p
     read(41, *) h1(j), J1(j, :)
     read(42, *) h3(j), J3(j, :)
  end do
  close(41)
  close(42)
  h2 = 0.0
  J2 = 0.0
  J2(1, 2) = 1.0
  J2(2, 1) = 1.0
  
  Z1 = 0.0
  Z2 = 0.0
  Z3 = 0.0
  do i = 0, 1
     x(1) = 2*i - 1
     do j = 0, 1
        x(2) = 2*j - 1
        do k = 0, 1
           x(3) = 2*k - 1
           do l = 0, 1
              x(4) = 2*l - 1
              call PartitionTerm(x, p, h1, J1, f)
              Z1 = Z1 + f
              call PartitionTerm(x, p, h2, J2, f)
              Z2 = Z2 + f
              call PartitionTerm(x, p, h3, J3, f)
              Z3 = Z3 + f
           end do
        end do
     end do
  end do
  
  write(*,*) Z1, Z2, Z3
  
  R = 1.0
  do l = 1, n
     R = R*Z1/Z3
     x = samples(l, :)
     if ( l > 10000 ) then
        call PartitionTerm(x, p, h2, J2, p1)
     else
        call PartitionTerm(x, p, h1, J1, p1)     
     end if
     call PartitionTerm(x, p, h3, J3, p2)
     R = R*p2/p1
     write(*,*), x
     write(*,*), l, R, p2, p1
     read(*, *)
  end do
  write(*,*) R
end program Ratiotest

  subroutine PartitionTerm(x, p, h, JJ, f)
    implicit none
    integer, intent(in) :: x(p)
    integer, intent(in) :: p
    real, intent(in) :: h(p)
    real, intent(in) :: JJ(p, p)
    real, intent(out) :: f
    integer i, j
    
    f = 0.0
    do i = 1, p
       f = f + x(i)*h(i)
       do j = 1, i-1
          f = f + x(i)*x(j)*JJ(i, j)
       end do
    end do
    !write(*,*), 'x', x
    !write(*,*), 'f', f   
    f = exp(f)
  end subroutine PartitionTerm
