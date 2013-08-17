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
    integer n_newborn
    real ls, us, ds
    integer Cindex
    character(30) filename

    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()
    open (unit = 17, file="randomnumber", action="write")

    t = 0.0
    tp = 0.0
    output_index = 0
    n_newborn = 0

    do while (t < tend)
       if (t .ge. tp) then
          call output_to_file(output_index)
          call output_to_file_newborn(output_index)
          output_index = output_index + 1
          tp = tp + tpinc
          print *, t
       end if

       do i = 1, NPool
          lsize = CellPool(i)%Csize
          lage = CellPool(i)%Cage
          lx(1) = CellPool(i)%mRNA
          lx(2) = CellPool(i)%nRibsome

          call getrate(lx, lage, a)
          do j = 1, lReac
             r(j) = poidev(a(j)*timestep)
             if (r(2) > lx(1) ) then
                r(2) = lx(1)
             end if
             if (r(4) > lx(2) ) then
                r(4) = lx(2)
             end if
             lx = lx + nu(:, j)*r(j)
             !print *, j, r(j), a(j)
          end do
          call checkx(lx, a, r, is_nag)          
          lsize = max(lsize, lx(2)/rho)
          lage = lage + timestep
          CellPool(i)%Csize = lsize
          CellPool(i)%Cage = lage
          CellPool(i)%mRNA = lx(1)
          CellPool(i)%nRibsome = lx(2)
          
          call check_mitosis(lsize, lage, m_flag)
          if ( m_flag .eq. 1 ) then
             call cell_division(i)
             ! Collect the newborn cells
             n_newborn = mod(n_newborn, NCollect) + 1
             NewbornPool(n_newborn) = CellPool(i)
          end if
       end do
       t = t + timestep

    end do
    close(unit=100)

  end subroutine cell_grow

end module cellgrow
