module setting
  integer, parameter :: L = 600
  integer, parameter :: D = 100
  integer, parameter :: b = 20
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 5000.0
  real p0(L)
  real v0(L)
  real p1(L)
  real, parameter :: v0m = 8.0
  real, parameter :: v_high = 1.0
  real, parameter :: vdiff_high = 0.2
  real, parameter :: vdiff_low = 0.2
  real TGFbeta(1-b:L+b)

contains
  subroutine output_to_file(cell_matrix, index)
    implicit none

    integer, intent(in) :: cell_matrix(0:L+1, D)
    integer, intent(in) :: index
    character(30) filename
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, L
       do j = 1, D
          write(11, '(I5)', advance="no"), cell_matrix(i,j)
       end do
       write(11, '(F15.2)', advance="no"), p0(i)
       write(11, '(F15.2)', advance="no"), TGFbeta(i)
       write(11, '(F15.2)', advance="no"), v0(i)
       write(11, *)
    end do
    close(11)
  end subroutine output_to_file

end module setting
