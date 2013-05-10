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
    real a(lReac), lx(lSpec)

    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    tm = 0.0
    output_index = 0

    do while (t < tend)
       
       !print *, t
       
       if (t .ge. tp) then
          call output_to_file(output_index)
          output_index = output_index + 1
          tp = tp + tpinc
          print *, t
       end if

       if (t .ge. tm) then
          do i = 1, NPool
             if ( CellPool(i)%Cage > tminc ) then
                CellPool(i)%Csize_old = CellPool(i)%Csize
             end if
          end do
          tm = tm + tminc
       end if

       do i = 1, NPool
          lsize = CellPool(i)%Csize
          lage = CellPool(i)%Cage
          lx(1) = CellPool(i)%mRNA
          lx(2) = CellPool(i)%nRibsome

          call getrate(lx, lage, a)
          do j = 1, lReac
             r(j) = poidev(a(j)*timestep)
             lx = lx + nu(:, j)*r(j)
             !print *, j, r(j), a(j)
          end do
          call checkx(lx, a, r, is_nag)          
          lsize = max(lsize, lx(2))
          lage = lage + timestep
          CellPool(i)%Csize = lsize
          CellPool(i)%Cage = lage
          CellPool(i)%mRNA = lx(1)
          CellPool(i)%nRibsome = lx(2)
          
          call check_mitosis(lsize, lage, m_flag)
          if ( m_flag .eq. 1 ) then
             call cell_division(i)
          end if
       end do
       t = t + timestep

       !print *, t, CellPool(1)%Csize, CellPool(1)%Cage, & 
       !     CellPool(1)%Crate, CellPool(1)%mRNA, CellPool(1)%nRibsome 
       !if (CellPool(1)%Cage.eq.0.0) read(*,*)

    end do
    close(unit=100)

  end subroutine cell_grow

end module cellgrow
