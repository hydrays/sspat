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
    integer n_mitosis, n_newborn, n_newborn2
    integer num_livecell, first_record
    real prodiv(NSize)
    real basediv1(NSize), basediv2(NSize)
    real ls, us, ds
    integer Cindex
    character(30) filename

    open (unit = 100, file='./out/logfile', action="write")
    call ran_seed(sequence=iseed)
    call init_cell_pool()
    open (unit = 17, file="randomnumber", action="write")


    t = 0.0
    tp = 0.0
    tm = 0.0
    output_index = 0
    n_mitosis = 0
    n_newborn = 0
    n_newborn2 = 0

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
                CellPool(i)%Csize_old1 = CellPool(i)%Csize_old2
                CellPool(i)%Csize_old2 = CellPool(i)%Csize
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
             if (r(2) > lx(1) ) then
                r(2) = lx(1)
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
             ! Collect the size of mitosis cell
             if ( t > tc .and. n_mitosis < NCollect) then
                n_mitosis = n_mitosis + 1
                MitosisPool(n_mitosis) = CellPool(i)
             end if
             call cell_division(i)
             ! Collect the newborn cells
             if ( t > tc .and. n_newborn < NCollect) then
                n_newborn = n_newborn + 1
                NewbornPool(n_newborn) = CellPool(i)
             end if

             if ( t > tc .and. n_newborn2 < NCollect2) then
                n_newborn2 = n_newborn2 + 1
                NewbornPool2(n_newborn2) = CellPool(i)
             end if

          end if
       end do
       t = t + timestep

       !print *, t, CellPool(1)%Csize, CellPool(1)%Cage, & 
       !     CellPool(1)%Crate, CellPool(1)%mRNA, CellPool(1)%nRibsome 
       !if (CellPool(1)%Cage.eq.0.0) read(*,*)

       ! Collect new born cell

          
    end do
    close(unit=100)

    ! Output mitosis cell to file
    if ( n_mitosis < NCollect ) then
       print *, "Not enough mitosis cells are collected..."
       read(*,*)
    else
       call output_to_file(-1)
    end if

    ! Output newborn cell to file
    if ( n_newborn < NCollect ) then
       print *, "Not enough newborn cells are collected..."
       read(*,*)
    else
       call output_to_file(-2)
    end if

  end subroutine cell_grow

end module cellgrow
