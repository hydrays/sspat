module setting
  integer, parameter :: L = 200
  integer, parameter :: H = 100
  integer, parameter :: b = 10
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 5000.0
  real p0(L)
  real p1(L)
  real, parameter :: v = 1.0
  real, parameter :: D = 1.0
  real TGFbeta(1-b:L+b)

contains
  subroutine output_to_file(cmat, index)
    implicit none

    integer, intent(in) :: cmat(0:L+1, H)
    integer, intent(in) :: index
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, L
       do j = 1, H
          write(11, '(I5)', advance="no"), cmat(i,j)
       end do
       write(11, '(F15.2)', advance="no"), p0(i)
       write(11, '(F15.2)', advance="no"), TGFbeta(i)
       write(11, *)
    end do
    close(11)
  end subroutine output_to_file

end module setting
