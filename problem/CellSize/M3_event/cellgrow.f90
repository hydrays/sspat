module cellgrow
  use setting
  implicit none

contains
  subroutine cell_grow()
    use random
    use setting
    implicit none

    real t, tp, tm, u
    real lsize, lage
    integer output_index, i, j, r(lReac)
    integer is_nag, m_flag
    real a(lReac), lx(lSpec), rr(lReac)
    integer n_mitosis, n_newborn, n_newborn2
    integer num_livecell, first_record
    real prodiv(NSize)
    real basediv1(NSize), basediv2(NSize)
    real yp(lSpec), f2(lSpec), f3(lSpec) 
    real f4(lSpec), f5(lSpec), f1(lSpec)
    real ls, us, ds, ch
    integer Cindex
    character(30) filename

    call ran_seed(sequence=iseed)
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    tm = 0.0
    output_index = 0
    n_mitosis = 0
    n_newborn = 0
    n_newborn2 = 0

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          call output_to_file_newborn(output_index)
          output_index = output_index + 1
          tp = tp + tpinc
          print *, t
       end if

       if (t .ge. tm) then
          do i = 1, NPool
             if ( CellPool(i)%Cage > tminc ) then
                CellPool(i)%Csize_old1 = CellPool(i)%Csize_old2
                CellPool(i)%Csize_old2 = CellPool(i)%Csize
             end if
          end do
          tm = tm + tminc
       end if

       do i = 1, NPool
          lage = CellPool(i)%Cage
          lx(1) = CellPool(i)%mRNA
          lx(2) = CellPool(i)%Csize

          ! ODE update
          ! call getrate(lx, lage, a)
          ! lx(1) = lx(1) + (a(1) - a(2))*timestep
          ! lx(2) = lx(2) + max(0.0, a(3) - a(4))*timestep

          ! RK4 update
          call getderivative(lage, lx, yp)
          ch = timestep/4.0
          f5 = lx + ch*yp
          call getderivative(lage+ch, f5, f1 )
          ch = 3.0*timestep/32.0
          f5 = lx + ch*( yp + 3.0*f1 )
          call getderivative(lage+3.0*timestep/8.0, f5, f2 )
          ch = timestep/2197.0
          f5 = lx + ch*( 1932.0*yp &
               + ( 7296.0*f2 - 7200.0*f1 ) )
          call getderivative (lage+12.0*timestep/13.0, f5, f3 )
          ch = timestep/4104.0
          f5 = lx + ch*( ( 8341.0*yp &
               - 845.0*f3 ) + ( 29440.0*f2 &
               - 32832.0*f1 ) )
          call getderivative ( lage+timestep, f5, f4 )
          ch = timestep/20520.0
          f1 = lx + ch*( (-6080.0*yp &
               + ( 9295.0*f3 - 5643.0*f4 ) ) &
               + ( 41040.0*f1 - 28352.0*f2 ) )
          call getderivative ( lage + timestep/2.0, f1, f5 )
          !  Ready to compute the approximate solution at T+H.
          ch = timestep / 7618050.0
          lx = lx + ch * ( ( 902880.0 * yp &
               + ( 3855735.0 * f3 - 1371249.0*f4 ) ) &
               + ( 3953664.0 * f2 + 277020.0*f5 ) )
          
          if ( lx(2) .le. lx(1) ) then
             CellPool(i)%event = CellPool(i)%event + timestep*lx(2)
          end if
          
          ! Stochastic update
          ! do j = 1, lReac
          !    !r(j) = poidev(a(j)*timestep)
          !    rr(j) = a(j)*timestep
          !    ! if (r(2) > lx(1) ) then
          !    !    r(2) = lx(1)
          !    ! end if
          !    ! if (r(4) > lx(2) ) then
          !    !    r(4) = lx(2)
          !    ! end if
          !    lx = lx + nu(:, j)*rr(j)
          !    !print *, j, r(j), a(j)
          ! end do
          ! call checkx(lx, a, r, is_nag)          
          lage = lage + timestep
          CellPool(i)%Cage = lage
          CellPool(i)%mRNA = lx(1)
          CellPool(i)%Csize = lx(2)

          call check_mitosis(i, m_flag)
          if ( m_flag .eq. 1 ) then
             ! Collect the size of mitosis cell
             call cell_division(i)
             ! Collect the newborn cells
             n_newborn = mod(n_newborn, NPool) + 1
             NewbornPool(n_newborn) = CellPool(i)
          end if
       end do
       t = t + timestep

       !print *, t, CellPool(1)%Csize, CellPool(1)%Cage, & 
       !     CellPool(1)%Crate, CellPool(1)%mRNA, CellPool(1)%nRibsome 
       !if (CellPool(1)%Cage.eq.0.0) read(*,*)

       ! Collect new born cell
    end do
  end subroutine cell_grow
end module cellgrow
