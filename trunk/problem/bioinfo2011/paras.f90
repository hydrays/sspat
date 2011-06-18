module paras  
  use nrtype

  implicit none
  integer(I4B), parameter :: Nsample=100
  integer(I4B), parameter :: N=1000
  real(kind=8), parameter :: mu0=1e-3
  integer(I4B), parameter :: m=4
  integer(I4B), parameter :: mplusone=m+1
  integer(I4B), parameter :: msminusone=mplusone*mplusone-1  
  integer(I4B), parameter :: ms=mplusone*mplusone
  real(kind=8), parameter :: delta_s = 0.01
  real(kind=8), parameter :: delta_mu = 0.0

  integer(I4B) nu(0:m, 0:msminusone)

  integer(I4B), parameter :: MSN=1
  real(kind=8), parameter :: null_a=50.0*MSN
  integer(I4B), parameter :: Nc=10

contains

  subroutine getrate(x, s, mu, a)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    real(kind=8), intent(out) :: a(0:msminusone)
    real(kind=8), intent(in) :: s(0:m)
    real(kind=8), intent(in) :: mu(0:m)
    integer(I4B) i, j
    real(kind=8) sx(0:m), sum_sx
    sx = s * x
    sum_sx = sum(sx)
    do i = 0, m
       do j = 0, m
          if (i .eq. j) then ! mutation rate for i
             a(i*(mplusone) + j) = mu(i) * x(i)
          else ! i-type die and j-type born
             a(i*(mplusone) + j) = x(i) * (sx(j)/sum_sx)
          end if
       end do

    end do
  end subroutine getrate

  subroutine makenu(nu)
    implicit none
    integer(I4B) nu(0:m, 0:msminusone)
    integer(I4B) i, j
    nu = 0
    do i = 0, m
       do j = 0, m
          if (i .eq. j) then ! mutation rate for i
             nu(i, i*(mplusone) + j) = -1
             nu(i+1, i*(mplusone) + j) = 1
          else ! i-type die and j-type born
             nu(i, i*(mplusone) + j) = -1
             nu(j, i*(mplusone) + j) = 1
          end if
       end do
    end do
  end subroutine makenu
  
  subroutine makes(s, ds)
    implicit none
    real(kind=8), intent(out) :: s(0:m)    
    real(kind=8), intent(in) :: ds    
    integer(I4B) j
    s(0) = 1.0
    do j=1,m
       s(j) = s(j-1)*(1.0+ds)
    end do
  end subroutine makes

  subroutine makemu(mu)
    implicit none
    real(kind=8), intent(out) :: mu(0:m)    
    integer(I4B) j
    mu(0) = mu0
    do j=1,m
       mu(j) = mu(0) + delta_mu*j
    end do
  end subroutine makemu

  subroutine makermax(x, rmax)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    integer(I4B), intent(out) :: rmax(0:msminusone)
    integer(I4B) i, j
    rmax = 0
    do i = 0, m
       do j = 0, m
          if (i .eq. j) then ! mutation rate for i
             rmax(i*(mplusone) + j) = x(i)
          else ! i-type die and j-type born
             rmax(i*(mplusone) + j) = x(i)
          end if
       end do
    end do
  end subroutine makermax

  subroutine checkx(x, is_nag)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    integer(I4B), intent(out) :: is_nag
    integer(I4B) i
    do i = 0, m
       if (x(i) < 0) then
          is_nag = 1
       end if
    end do
  end subroutine checkx

  subroutine partition(x, critical_index)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    integer(I4B), intent(out) :: critical_index(0:msminusone)
    integer(I4B) i,j
    critical_index = 1
    do i = 0, m
       do j = 0, m
          if ( (x(i)>Nc) .and. (x(j)>Nc) .and. (i .ne. j) ) then
             critical_index(i*(mplusone) + j) = 0
          end if
       end do
    end do
  end subroutine partition

end module paras
