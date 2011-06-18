module paras  
  use nrtype

  implicit none
  integer(I4B), parameter :: Nsample=100000
  integer(I4B), parameter :: N=10000
  real(kind=8), parameter :: mu0=1e-8
  integer(I4B), parameter :: m=1
  real(kind=8), parameter :: delta_s = 0.01
  real(kind=8), parameter :: delta_mu = 0.0
  real(kind=8), parameter :: delta_t_coarse = 0.02
  integer(I4B), parameter :: Nc=10

contains

  subroutine getrate(x, s, mu, a)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    real(kind=8), intent(out) :: a(0:m, 0:m)
    real(kind=8), intent(in) :: s(0:m)
    real(kind=8), intent(in) :: mu(0:m)
    integer(I4B) i, j
    real(kind=8) sx(0:m), sum_sx
    sx = s * x
    sum_sx = sum(sx)
    do i = 0, m
       do j = 0, m
          if (i .eq. j) then ! mutation rate for i
             a(i, j) = mu(i) * x(i)
          else ! i-type die and j-type born
             a(i, j) = x(i) * (sx(j)/sum_sx)
          end if
       end do

    end do
  end subroutine getrate

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

  subroutine partition(x, cmask)
    implicit none
    integer(I4B), intent(in) :: x(0:m)
    logical, intent(out) :: cmask(0:m, 0:m)
    integer(I4B) i,j
    cmask = .true.
    do i = 0, m
       do j = 0, m
          if ( (x(i)>Nc) .and. (x(j)>Nc) .and. (i .ne. j) ) then
             cmask(i, j) = .false.
          end if
       end do
    end do
  end subroutine partition

end module paras
