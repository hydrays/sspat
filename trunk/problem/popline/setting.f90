module setting
  integer, parameter :: L = 800
  integer, parameter :: H = 200
  real, parameter :: b = 8.0
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 50000.0
  real, parameter :: p1 = 0.3
  real, parameter :: v = 1.0
  real, parameter :: D = 1.0
  real p0
  real TGFbeta(1-b:L+b)
  
  type cell
     integer type
     real gene1
  end type cell

contains
  subroutine output_to_file(cmat, index)
    implicit none

    type(cell), intent(in) :: cmat(0:L+1, H)
    integer, intent(in) :: index
    character(30) filename, filename2
    integer i, j

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    WRITE(filename2,'(A7,I5.5,A4)') './out/g', index, '.dat'
    open (unit = 11, file=filename, action="write")
    open (unit = 12, file=filename2, action="write")

    do i = 1, L
       do j = 1, H
          write(11, '(I5)', advance="no"), cmat(i,j)%type
       end do
       write(11, '(F15.2)', advance="no"), TGFbeta(i)
       write(11, *)
    end do
    do i = 1, L
       do j = 1, H
          if (cmat(i,j)%type.eq.1) then
             write(12, '(I10, (F15.5))'), i, cmat(i,j)%gene1
          end if
       end do
    end do
    close(11)
    close(12)
  end subroutine output_to_file

  subroutine cell_event(cmat, i, j)
    use random
    implicit none
    type(cell), intent(inout) :: cmat(0:L+1, H)
    integer, intent(in) :: i, j
    integer k
    real u, p0
    if ( cmat(i,j)%type .eq. 1 ) then
       p0 = cmat(i,j)%gene1 + (1.0 - 2*cmat(i,j)%gene1) &
            / (1.0 + 0.01*TGFbeta(i))
       ! division
       do k=H, j+2, -1
          cmat(i, k) = cmat(i, k-1)
       end do
       call ran2(u)
       if ( u < p0 ) then
          ! SC -> 2SC
          cmat(i,j+1) = cmat(i,j)
       else
          ! SC -> 2TAC
          cmat(i,j)%type = 2
          cmat(i,j+1) = cmat(i,j) 
       end if
    else if ( cmat(i,j)%type .eq. 2 ) then
       ! division
       do k=H, j+2, -1
          cmat(i, k) = cmat(i, k-1)
       end do
       call ran2(u)
       if ( u < p1 ) then
          ! TAC -> 2TAC
          cmat(i, j+1) = cmat(i,j)
       else
          ! SC -> 2TAC
          cmat(i, j)%type = 3
          cmat(i, j+1) = cmat(i,j)
       end if
    else if ( cmat(i,j)%type .eq. 3 ) then
       ! death
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
    else
       ! do nothing
    end if
  end subroutine cell_event

  subroutine getTGFbeta(cmat)
    implicit none
    type(cell), intent(inout) :: cmat(0:L+1, H)
    integer i, j, k
    real temp_num

    TGFbeta = 0.
    do i = 1, L
       temp_num = 0.
       do j = 1, H
          if ( cmat(i,j)%type .eq. 3 ) then
             temp_num = temp_num + 1.
          end if
       end do
       if (temp_num > 0) then
        do k = 0, b-1
           if ( k .eq. 0) then
              TGFbeta(i) = TGFbeta(i) + temp_num
           else
              TGFbeta(i+k) = TGFbeta(i+k) + temp_num*exp(-real(k)/b)
              TGFbeta(i-k) = TGFbeta(i-k) + temp_num*exp(-real(k)/b)
           end if
        end do
     end if
  end do
  do k = 1, b
     TGFbeta(k) = TGFbeta(k) + TGFbeta(L+k)
     TGFbeta(L-k+1) = TGFbeta(L-k+1) + TGFbeta(1-k)
  end do

  end subroutine getTGFbeta
end module setting