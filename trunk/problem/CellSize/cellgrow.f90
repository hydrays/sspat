module cellgrow
  use setting
  implicit none

contains
  subroutine cell_grow()
    use random
    use setting
    implicit none

    real t, tp, u
    real lsize, lage
    integer output_index, i, j, r(lReac)
    integer is_nag, m_flag
    real a(lReac), lx(lSpec)

    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()

    t = 0.0
    tp = 0.0
    output_index = 0

    do while (t < tend)
       
       !print *, t
       
       if (t .ge. tp) then
          call output_to_file(output_index)
          output_index = output_index + 1
          tp = tp + tpinc
          print *, t
       end if

       do i = 1, NPool
          lsize = CellPool(i)%Csize
          lage = CellPool(i)%Cage
          lx(1) = CellPool(i)%mRNA
          lx(3) = CellPool(i)%nRibsome
          lx(2) = min(lx(1), lx(3))

          call getrate(lx, lage, a)
          do j = 1, lReac
             r(j) = poidev(a(j)*timestep)
             lx = lx + nu(:, j)*r(j)


             ! if ( j .eq. 3 ) then
             !    if ( r(j) > lx(1) ) then
             !       r(j) = lx(1)
             !    end if
             !    if ( r(j) + lx(2) > lx(3) ) then
             !       r(j) = lx(3) - lx(2)
             !    end if
             ! end if
             ! if ( j .eq. 4 .and. r(j) > lx(2) ) then
             !    r(j) = lx(2)
             ! end if

             !print *, j, r(j), a(j)
          end do
          call checkx(lx, is_nag)          
          if ( lsize < lx(3) ) then
             CellPool(i)%Crate = 0.0
          else
             CellPool(i)%Crate = (a(5) - a(6))
          end if
          lsize = max(lsize, lx(3))
          lage = lage + timestep
          CellPool(i)%Csize = lsize
          CellPool(i)%Cage = lage
          CellPool(i)%mRNA = lx(1)
          CellPool(i)%nActiveRibsome = lx(2)
          CellPool(i)%nRibsome = lx(3)
          
          call check_mitosis(lsize, lage, m_flag)
          if ( m_flag .eq. 1 ) then
             call cell_division(i)
          end if
       end do
       t = t + timestep

       print *, t, CellPool(1)%Csize, CellPool(1)%Cage, & 
            CellPool(1)%Crate, CellPool(1)%mRNA, &
            CellPool(1)%nActiveRibsome, CellPool(1)%nRibsome 
       if (CellPool(1)%Cage.eq.0.0) read(*,*)

    end do
    close(unit=100)

  end subroutine cell_grow

end module cellgrow
