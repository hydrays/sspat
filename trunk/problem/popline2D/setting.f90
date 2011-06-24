module setting
  integer, parameter :: L = 200
  integer, parameter :: H = 50
  real, parameter :: b = 8.0
  real, parameter :: delta_t = 0.001
  real, parameter :: tend = 50000.0
  real, parameter :: p1 = 0.3
  real, parameter :: v = 1.0
  real, parameter :: D = 1.0
  real p0
  real TGFbeta(1-b:L+b, 1-b:L+b)
  
  type cell
     integer type
     real gene1
  end type cell

contains
  subroutine output_to_file(cmat, index)
    implicit none

    type(cell), intent(in) :: cmat(0:L+1, 0:L+1, H)
    integer, intent(in) :: index
    character(30) filename, filename2
    integer ix, iy, j
  integer num_sc, num_tac, num_tdc, num_mc

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    WRITE(filename2,'(A7,I5.5,A4)') './out/g', index, '.dat'
    open (unit = 11, file=filename, action="write")
    open (unit = 12, file=filename2, action="write")

    do ix = 1, L
       do iy = 1, L
          num_sc = 0
          num_tac = 0
          num_tdc = 0
          num_mc = 0
          do j = 1, H
             if (cmat(ix,iy,j)%type.eq.1) then
                num_sc = num_sc + 1
             elseif (cmat(ix,iy,j)%type.eq.2) then
                num_tac = num_tac + 1
             elseif (cmat(ix,iy,j)%type.eq.3) then
                num_tdc = num_tdc + 1
             elseif (cmat(ix,iy,j)%type.eq.4) then
                num_mc = num_mc + 1
             end if
          end do
          write(11, '(7I12)'), ix, iy, num_sc, num_tac, &
               num_tdc, num_mc, num_sc + num_tac + num_tdc + num_mc
       end do
    end do
    do ix = 1, L
       do iy = 1, L
          do j = 1, H
             if (cmat(ix,iy,j)%type.eq.1) then
                write(12, '(2I10, (F15.5))'), ix, iy, cmat(ix,iy,j)%gene1
             end if
          end do
       end do
    end do
    close(11)
    close(12)
  end subroutine output_to_file

  subroutine cell_event(cmat, ix, iy, j)
    use random
    implicit none
    type(cell), intent(inout) :: cmat(0:L+1, 0:L+1, H)
    integer, intent(in) :: ix, iy, j
    integer k
    real u, p0
    if ( cmat(ix,iy,j)%type .eq. 1 ) then
       p0 = cmat(ix,iy,j)%gene1 + (1.0 - 2*cmat(ix,iy,j)%gene1) &
            / (1.0 + 0.01*TGFbeta(ix, iy))
       ! division
       do k=H, j+2, -1
          cmat(ix,iy,k) = cmat(ix,iy,k-1)
       end do
       call ran2(u)
       if ( u < p0 ) then
          ! SC -> 2SC
          cmat(ix,iy,j+1) = cmat(ix,iy,j)
       else
          ! SC -> 2TAC
          cmat(ix,iy,j)%type = 2
          cmat(ix,iy,j+1) = cmat(ix,iy,j) 
       end if
    else if ( cmat(ix,iy,j)%type .eq. 2 ) then
       ! division
       do k=H, j+2, -1
          cmat(ix,iy,k) = cmat(ix,iy,k-1)
       end do
       call ran2(u)
       if ( u < p1 ) then
          ! TAC -> 2TAC
          cmat(ix,iy,j+1) = cmat(ix,iy,j)
       else
          ! SC -> 2TAC
          cmat(ix,iy,j)%type = 3
          cmat(ix,iy,j+1) = cmat(ix,iy,j)
       end if
    else if ( cmat(ix,iy,j)%type .eq. 3 ) then
       ! death
       do k=j, H-1
          cmat(ix,iy,k) = cmat(ix,iy,k+1)
       end do
    else
       ! do nothing
    end if
  end subroutine cell_event

  subroutine getTGFbeta(cmat)
    implicit none
    type(cell), intent(inout) :: cmat(0:L+1,0:L+1,H)
    integer ix,iy, j, k, kx, ky
    real temp_num, dist

    TGFbeta = 0.
    do ix = 1, L
       do iy = 1, L
          temp_num = 0.
          do j = 1, H
             if ( cmat(ix,iy,j)%type .eq. 3 ) then
                temp_num = temp_num + 1.
             end if
          end do
          if (temp_num > 0) then
             do kx = -b, b
                do ky = -b, b
                   dist = sqrt(real(kx*kx + ky*ky))
                   TGFbeta(ix+kx,iy+ky) = &
                        TGFbeta(ix+kx,iy+ky) + temp_num*exp(-dist/real(b))
                end do
             end do
          end if
       end do
    end do
    do k = 1, b
       TGFbeta(1:L,k) = TGFbeta(1:L,k) + TGFbeta(1:L,L+k)
       TGFbeta(1:L,L-k+1) = TGFbeta(1:L,L-k+1) + TGFbeta(1:L,1-k)
       TGFbeta(k,1:L) = TGFbeta(k,1:L) + TGFbeta(L+k,1:L)
       TGFbeta(L-k+1,1:L) = TGFbeta(L-k+1,1:L) + TGFbeta(1-k,1:L)
    end do
  end subroutine getTGFbeta

end module setting
