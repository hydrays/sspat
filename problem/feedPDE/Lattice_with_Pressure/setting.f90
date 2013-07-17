module setting
  integer :: L
  integer :: H
  integer :: brange, iseed
  real :: tend, p1, v, difv, mutv, tm
  real :: fdgain1, scstick, prelax, tpinc
  real :: bd10
  integer :: useomp, is64bit
  real :: timestep
  integer :: npar

  type cell
     integer type
     real gene1
     real gene2
     real gene3
     real gene4
  end type cell
  type(cell), allocatable :: cmat(:,:)
  real, allocatable :: a(:)
  real, allocatable :: NT(:)
  real, allocatable :: NP(:)

  ! Nutrition that effect division rate
  real, allocatable :: Nutri(:), Nutri_old(:)
  real NutriGrowthRate
  real NutriConsumeRate
  real NutriDecayRate
  real NutriMobility
  real NutriTimestep
  real NutriKillrate

  integer, allocatable :: npack(:)
  integer, allocatable :: TDC(:)
  integer, allocatable :: SC(:)
  integer, allocatable :: TAC(:)
  integer, allocatable :: MC(:)

  namelist /xdata/ L, H, brange, tend, p1, v, difv, mutv, &
       fdgain1, scstick, prelax, iseed, tpinc, tm

  namelist /xdataomp/ useomp, is64bit, timestep, npar

  namelist /ndata/ NutriGrowthRate, NutriConsumeRate, &
       NutriDecayRate, NutriMobility, NutriTimestep, NutriKillrate

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)
    read(8, nml=xdataomp)
    read(8, nml=ndata)

    bd10 = 10.0/real(brange)

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'L = ', L
    write(*, '(a20, i10)'), 'H = ', H
    write(*, '(a20, i10)'), 'brange = ', brange
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, f10.2)'), 'tend = ', tend
    write(*, '(a20, f10.2)'), 'p1 = ', p1
    write(*, '(a20, f10.2)'), 'v = ', v
    write(*, '(a20, f10.2)'), 'difv = ', difv
    write(*, '(a20, f10.2)'), 'mutv = ', mutv
    write(*, '(a20, f10.2)'), 'fdgain1 = ', fdgain1
    write(*, '(a20, f10.2)'), 'scstick = ', scstick
    write(*, '(a20, f10.2)'), 'prelax = ', prelax
    write(*, '(a20, f10.2)'), 'tpinc = ', tpinc
    write(*, '(a20, f10.2)'), 'tm = ', tm
    write(*, '(a20, f10.2)'), 'NutriGrowthRate = ', NutriGrowthRate
    write(*, '(a20, f10.2)'), 'NutriConsumeRate = ', NutriConsumeRate
    write(*, '(a20, f10.2)'), 'NutriDecayRate = ', NutriDecayRate
    write(*, '(a20, f10.2)'), 'NutriMobility = ', NutriMobility
    write(*, '(a20, f10.2)'), 'NutriTimestep = ', NutriTimestep
    write(*, '(a20, f10.2)'), 'NutriKillrate = ', NutriKillrate

    if (useomp.eq.1) then
       write(*, '(a)'), 'OpenMP parallel in use!'
       write(*, '(a20, I10)'), 'is64bit = ', is64bit
       write(*, '(a20, f10.2)'), 'timestep = ', timestep
       write(*, '(a20, I10)'), 'npar = ', npar
    end if

    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'L,', L
    write(9, '(a20, i10)'), 'H,', H
    write(9, '(a20, i10)'), 'brange,', brange
    write(9, '(a20, i10)'), 'iseed,', iseed
    write(9, '(a20, f10.2)'), 'tend,', tend
    write(9, '(a20, f10.2)'), 'p1,', p1
    write(9, '(a20, f10.2)'), 'v,', v
    write(9, '(a20, f10.2)'), 'difv,', difv
    write(9, '(a20, f10.2)'), 'mutv,', mutv
    write(9, '(a20, f10.2)'), 'fdgain1,', fdgain1
    write(9, '(a20, f10.2)'), 'scstick,', scstick
    write(9, '(a20, f10.2)'), 'prelax,', prelax
    write(9, '(a20, f10.2)'), 'tpinc,', tpinc
    write(9, '(a20, f10.2)'), 'tm,', tm
    write(9, '(a20, i10)'), 'useomp,', useomp
    write(9, '(a20, i10)'), 'is64bit,', is64bit
    write(9, '(a20, f10.2)'), 'timestep,', timestep
    write(9, '(a20, I10)'), 'npar,', npar
    write(9, '(a20, f10.2)'), 'NutriGrowthRate,', NutriGrowthRate
    write(9, '(a20, f10.2)'), 'NutriConsumeRate,', NutriConsumeRate
    write(9, '(a20, f10.2)'), 'NutriDecayRate,', NutriDecayRate
    write(9, '(a20, f10.2)'), 'NutriMobility,', NutriMobility
    write(9, '(a20, f10.2)'), 'NutriTimestep,', NutriTimestep
    write(9, '(a20, f10.2)'), 'NutriKillrate,', NutriKillrate
    close(8)
    close(9)

  end subroutine read_xdata

  subroutine init_cell_pool()
    use random
    implicit none
    integer temp_num
    real u
    integer i, j

    allocate(cmat(0:L+1,H))
    allocate(a(1:L))
    allocate(NT(1:L))
    allocate(NP(1:L))
    allocate(npack(0:L+1))
    allocate(TDC(0:L+1))
    allocate(SC(0:L+1))
    allocate(MC(0:L+1))
    allocate(TAC(0:L+1))
    allocate(Nutri(0:L+1))
    allocate(Nutri_old(0:L+1))

!    cmat(200, 1)%type = 1
    cmat(1:L, 1)%type = 1
    cmat(1:L, 2:3)%type = 2
    cmat(1:L, 3:10)%type = 3
    do i = 1, L
       call ran2(u)
       cmat(i, 1)%gene1 = scstick
       cmat(i, 1)%gene2 = prelax
       cmat(i, 1)%gene3 = fdgain1
       cmat(i, 1)%gene4 = 10!+2.0*(u-0.5)
    end do
    cmat(0, :) = cmat(L, :)
    cmat(L+1, :) = cmat(1, :)

    npack = 0
    TDC = 0
    SC = 0
    MC = 0
    TAC = 0
    do i = 1, L
       do j = 1, H
          if ( cmat(i,j)%type .ne. 0 ) then
             npack(i) = npack(i) + 1
          end if
          if ( cmat(i,j)%type .eq. 1 ) then
             SC(i) = SC(i) + 1
          else if ( cmat(i,j)%type .eq. 2 ) then
             TAC(i) = TAC(i) + 1
          else if ( cmat(i,j)%type .eq. 3 ) then
             TDC(i) = TDC(i) + 1
          else if ( cmat(i,j)%type .eq. 4 ) then
             MC(i) = MC(i) + 1
          end if
       end do
    end do
    npack(0) = npack(L)
    npack(L+1) = npack(1)
    TDC(0) = TDC(L)
    TDC(L+1) = TDC(1)
    SC(0) = SC(L)
    SC(L+1) = SC(1)
    TAC(0) = TAC(L)
    TAC(L+1) = TAC(1)
    MC(0) = MC(L)
    MC(L+1) = MC(1)

    NP = 0.0
    NT = 0.0
    do i = 1, L
       call random_number(u)
       u = -log(u)
       call Update_Rate(i)
       NP(i) = u
    end do

    ! Nutrition inital is high
    Nutri(0:L+1) = 1.0

  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename, filename2
    integer i, j
    
    integer k, shift_i
    real TGFbeta, p0

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    WRITE(filename2,'(A7,I5.5,A4)') './out/g', index, '.dat'
    open (unit = 11, file=filename, action="write")
    open (unit = 12, file=filename2, action="write")

    do i = 1, L
       do j = 1, H
          write(11, '(I5)', advance="no"), cmat(i,j)%type
       end do
       write(11, '(f10.2)', advance="no"), Nutri(i)
       write(11, '(I5)', advance="no"), SC(i)
       write(11, '(I5)', advance="no"), TAC(i)
       write(11, '(I5)', advance="no"), TDC(i)
       write(11, '(I5)', advance="no"), MC(i)
       write(11, '(I5)', advance="no"), npack(i)

       TGFbeta = 0.0
       do k = -brange, brange
          shift_i = k + i
          if ( shift_i .le. 0 ) then
             shift_i = shift_i + L
          else if ( shift_i > L ) then
             shift_i = shift_i - L
          end if
          TGFbeta = TGFbeta + &
               bd10*TDC(shift_i)*exp(-real(abs(k))/brange)
       end do
       p0 = 0.2 + 0.6 / (1.0 + 0.01*TGFbeta)

       write(11, '(f10.2)', advance="no"), p0
       write(11, '(f10.2)', advance="no"), ((25*Nutri(i))**2)/(1.0+((25*Nutri(i))**2))
       write(11, *)
    end do
    do i = 1, L
       do j = 1, H
          if (cmat(i,j)%type.eq.1) then
             write(12, '(I10, 4(F15.5))'), i, &
                  cmat(i,j)%gene1, cmat(i,j)%gene2, &
                  cmat(i,j)%gene3, cmat(i,j)%gene4
          end if
       end do
    end do
    close(11)
    close(12)
  end subroutine output_to_file

  subroutine cell_event(i)
    use random
    implicit none
    integer, intent(in) :: i
    integer j, k, m, shift_i
    real u, p0, TGFbeta, Pa, ENutri
    real vr, vl
    type(cell) new_cell
    real u1, ua, ud

    !Effective Nutrition
    !if ( Nutri(i) > 1.0 ) then
    !ENutri = min(5.0, Nutri(i))
    !end if
    !if ( Nutri(i) < 1.0 ) then
    !   ENutri = 0.0
    !end if
    !ENutri = ENutri / 5.0
    ! Pa is the probability of accepting the division of 
    ! stem cell and mutation cell effected by the nutrition
    ! supply from the basal layer, Nutri.
    ! Pa controls the division rate of SC and MC.
    !Pa = 2.0*ENutri/(1.0+ENutri)
    Pa = ((25*Nutri(i))**2)/(1.0+((25*Nutri(i))**2))

    call ran2(u)
    u = u*a(i)
    !if (npack(i).ne.0) then
       vr = max(0.0, difv*real(npack(i)-npack(i+1)))
       vl = max(0.0, difv*real(npack(i)-npack(i-1)))
       !vr = difv*real(npack(i))
       !vl = difv*real(npack(i))
    !else
    !   vr = 0.0
    !   vl = 0.0
    !end if
!!$    print *, 'event happen at ', i
!!$    print *, 'u', u
!!$    print *, 'a(i)', a(i)
!!$    print *, 'npack(i)', npack(i-1:i+1)
!!$    print *, 'vr, vl', vr, vl

    u = u - vr
    if ( u < 0 ) then
       call ran2(u1)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 4', j, u1
          read(*,*)
       end if

       if ( cmat(i, j)%type .eq. 1 .or. cmat(i,j)%type .eq. 4) then
          call ran2(u1)
          if ( u1 < cmat(i,j)%gene1 ) then
             return
          end if
       end if

       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i + 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1

       if (new_cell%type .eq. 1) then
          SC(i) = SC(i) - 1
          SC(m) = SC(m) + 1
       else if (new_cell%type .eq. 2) then
          TAC(i) = TAC(i) - 1
          TAC(m) = TAC(m) + 1
       else if (new_cell%type .eq. 3) then
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       else if (new_cell%type .eq. 4) then
          MC(i) = MC(i) - 1
          MC(m) = MC(m) + 1
       end if
       return
    end if
    u = u - vl
    if ( u < 0 ) then
       call ran2(u1)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 5', j, u1, npack(i)
          read(*,*)
       end if
       !print *, 'move left at height j', i, j
       if ( cmat(i, j)%type .eq. 1 .or. cmat(i,j)%type .eq. 4) then
          call ran2(u1)
          if ( u1 < cmat(i,j)%gene1 ) then
             return
          end if
       end if
       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i - 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1
       if (new_cell%type .eq. 1) then
          SC(i) = SC(i) - 1
          SC(m) = SC(m) + 1
       else if (new_cell%type .eq. 2) then
          TAC(i) = TAC(i) - 1
          TAC(m) = TAC(m) + 1
       else if (new_cell%type .eq. 3) then
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       else if (new_cell%type .eq. 4) then
          MC(i) = MC(i) - 1
          MC(m) = MC(m) + 1
       end if
       return
    end if
    do j = 1, npack(i)
       u = u - v
       if ( u < 0 ) then
          if ( cmat(i,j)%type .eq. 1 ) then
             call ran2(ua)
             if (ua .le. Pa) then
                TGFbeta = 0.0
                do k = -brange, brange
                   shift_i = k + i
                   if ( shift_i .le. 0 ) then
                      shift_i = shift_i + L
                   else if ( shift_i > L ) then
                      shift_i = shift_i - L
                   end if
                   TGFbeta = TGFbeta + &
                        bd10*TDC(shift_i)*exp(-real(abs(k))/brange)
                end do
                p0 = cmat(i,j)%gene2 + (1.0 - 2.0*cmat(i,j)%gene2) &
                     / (1.0 + cmat(i,j)%gene3*TGFbeta)
                !p0 = 0.2 + 0.6 / (1.0 + 0.01*TGFbeta)
                !p0 = p0*(1.0-real(j)/40.0)
                !print *, 'p0', p0
                ! division
                do k=H, j+2, -1
                   cmat(i, k) = cmat(i, k-1)
                end do
                call ran2(u1)
                if ( u1 < p0 ) then
                   ! SC -> 2SC
                   cmat(i,j+1) = cmat(i,j)
                   SC(i) = SC(i) + 1
                else
                   ! SC -> 2TAC
                   cmat(i,j)%type = 2
                   cmat(i,j+1) = cmat(i,j)
                   SC(i) = SC(i) - 1
                   TAC(i) = TAC(i) + 2
                end if
                npack(i) = npack(i) + 1
             else
                !call ran2(ud)
                !if ( Nutri(i) < 0.001 ) then
                !   ud = ud * 0.02
                !end if
                !if ( ud < NutriKillrate ) then
                !if ( Pa .eq. 0.0 ) then
                if ( Nutri(i) < 0.05*0.04 ) then
                   ! death
                   do k=j, H-1
                      cmat(i, k) = cmat(i, k+1)
                   end do
                   SC(i) = SC(i) - 1
                   npack(i) = npack(i) - 1
                end if
             end if
          else if ( cmat(i,j)%type .eq. 2 ) then
             ! division
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             call ran2(u1)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i, j+1) = cmat(i,j)
                TAC(i) = TAC(i) + 1
             else
                ! TAC -> 2TDC
                cmat(i, j)%type = 3
                cmat(i, j+1) = cmat(i,j)
                TAC(i) = TAC(i) - 1
                TDC(i) = TDC(i) + 2
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 3 ) then
             ! death
             do k=j, H-1
                cmat(i, k) = cmat(i, k+1)
             end do
             TDC(i) = TDC(i) - 1
             npack(i) = npack(i) - 1
          else if ( cmat(i,j)%type .eq. 4 ) then
             call ran2(ua)
             if (ua .le. Pa) then
                do k=H, j+2, -1
                   cmat(i, k) = cmat(i, k-1)
                end do
                ! MC -> 2MC
                cmat(i,j+1) = cmat(i,j)
                MC(i) = MC(i) + 1
                npack(i) = npack(i) + 1
             else
                !call ran2(ud)
                !if ( Nutri(i) < 0.001 ) then
                !   ud = ud * 0.02
                !end if
                !if ( ud < NutriKillrate ) then
                !if ( Pa .eq. 0.0 ) then
                if ( Nutri(i) < 0.05*0.04 ) then
                   ! death
                   do k=j, H-1
                      cmat(i, k) = cmat(i, k+1)
                   end do
                   MC(i) = MC(i) - 1
                   npack(i) = npack(i) - 1
                end if
             end if
          else
             ! do nothing
          end if
          return
       end if
    end do
    if ( mutv .ne. 0.0 ) then
       do j = 1, npack(i)
          u = u - mutv
          if ( u < 0 ) then
             if ( cmat(i,j)%type .eq. 1 ) then
                ! mutation
                call ran2(u1)
                cmat(i,j)%gene1 = cmat(i,j)%gene1 + (u1-0.5)*0.1
                if ( cmat(i,j)%gene1 > 0.995) then
                   cmat(i,j)%gene1 = 0.995
                else if ( cmat(i,j)%gene1 < 0.0) then
                   cmat(i,j)%gene1 = 0.0
                end if
!!$                u1 = par_uni(kpar)
!!$                cmat(i,j)%gene2 = cmat(i,j)%gene2 + (u1-0.5)*0.01
!!$                if ( cmat(i,j)%gene2 > 0.5) then
!!$                   cmat(i,j)%gene2 = 1.0
!!$                else if ( cmat(i,j)%gene2 < 0.0) then
!!$                   cmat(i,j)%gene2 = 0.0
!!$                end if
!!$                u1 = par_uni(kpar)
!!$                cmat(i,j)%gene3 = cmat(i,j)%gene3 + (u1-0.5)*0.0001
!!$                if ( cmat(i,j)%gene2 > 0.5) then
!!$                   cmat(i,j)%gene2 = 1.0
!!$                else if ( cmat(i,j)%gene2 < 0.0) then
!!$                   cmat(i,j)%gene2 = 0.0
!!$                end if
             end if
             return
          end if
       end do
    end if
    print *, "not suppose to be here!"
    stop
  end subroutine cell_event

  subroutine cell_restack(i)
    implicit none
    integer, intent(in) :: i
    integer j
    type(cell) temp
    do j = npack(i), 2, -1
       if ( cmat(i,j)%type .eq. 1 ) then
          temp = cmat(i,j-1)
          cmat(i,j-1) = cmat(i,j)
          cmat(i,j) = temp
       end if
    end do
  end subroutine cell_restack

  subroutine cell_stat(t)
    implicit none

    real, intent(in) :: t
    integer i, j, k
    integer num_sc, num_tac, num_tdc, num_mc

    num_sc = 0
    num_tac = 0
    num_tdc = 0
    num_mc = 0
    do i = 1, L
       do j = 1, H
          if (cmat(i,j)%type.eq.1) then
             num_sc = num_sc + 1
          elseif (cmat(i,j)%type.eq.2) then
             num_tac = num_tac + 1
          elseif (cmat(i,j)%type.eq.3) then
             num_tdc = num_tdc + 1
          elseif (cmat(i,j)%type.eq.4) then
             num_mc = num_mc + 1
          end if
       end do
    end do

    write(*, '(5(F10.2))'), t, real(num_sc)/L, &
         real(num_tac)/L, real(num_tdc)/L, real(num_mc)/L
!!$  do i = -b, L+1+b
!!$     write(*, '(F8.2)', advance="no"), TGFbeta(i)
!!$  end do
!!$  write(*, *)
    write(100, '(5(F10.2))'), t, real(num_sc)/L, &
         real(num_tac)/L, real(num_tdc)/L, real(num_mc)/L
  end subroutine cell_stat

  subroutine Perodic_BC(k)
    implicit none
    integer, intent(in) :: k

    if ( k .eq. 1 ) then
       cmat(L, :) = cmat(0, :)
       cmat(L+1, :) = cmat(1, :)
       npack(L) = npack(0)
       npack(L+1) = npack(1)
       TDC(L) = TDC(0)
       TDC(L+1) = TDC(1)
       SC(L) = SC(0)
       SC(L+1) = SC(1)
       TAC(L) = TAC(0)
       TAC(L+1) = TAC(1)
       MC(L) = MC(0)
       MC(L+1) = MC(1)
    end if
    if ( k .eq. 2 ) then
       cmat(L+1, :) = cmat(1, :)
       npack(L+1) = npack(1)
       TDC(L+1) = TDC(1)
       SC(L+1) = SC(1)
       TAC(L+1) = TAC(1)
       MC(L+1) = MC(1)
    end if
    if ( k .eq. L ) then
       cmat(1, :) = cmat(L+1, :)
       cmat(0, :) = cmat(L, :)
       npack(1) = npack(L+1)
       npack(0) = npack(L)
       TDC(1) = TDC(L+1)
       TDC(0) = TDC(L)
       SC(1) = SC(L+1)
       SC(0) = SC(L)
       TAC(1) = TAC(L+1)
       TAC(0) = TAC(L)
       MC(1) = MC(L+1)
       MC(0) = MC(L)
    end if
    if ( k .eq. L-1 ) then
       cmat(0, :) = cmat(L, :)
       npack(0) = npack(L)
       TDC(0) = TDC(L)
       SC(0) = SC(L)
       TAC(0) = TAC(L)
       MC(0) = MC(L)
    end if
  end subroutine Perodic_BC

  subroutine Next_Reaction(k, tau)
    implicit none
    integer, intent(out) :: k
    real, intent(out) :: tau
    real tau_temp
    integer i

    tau = huge(0.0)
    k = 0
    do i = 1, L
       if ( a(i) > 0.0 ) then
          tau_temp = ( NP(i) - NT(i) ) / a(i)
          if ( tau_temp < tau) then
             tau = tau_temp
             k = i
          end if
       else
          !          print *, 'no reaction -------------- ', i, k, a(i)
          !          read(*,*)
       end if
    end do
    !    print *, 'k', k
    !    print *, 'a', a
    if ( k .le. 0 ) then
       write(*,*), 'error'
       read(*,*)
    end if
    if ( tau < 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, 'NP-NT', NP - NT
       read(*,*)
    end if
  end subroutine Next_Reaction

  subroutine Update_Rate(i)
    implicit none
    integer, intent(in) :: i
    real vr, vl
    vr = max(0.0, difv*real(npack(i)-npack(i+1)))
    vl = max(0.0, difv*real(npack(i)-npack(i-1)))
    !vr = difv*real(npack(i))
    !vl = difv*real(npack(i))
    if (npack(i).eq.0) then
       !       print *, vr, vl
    end if
    a(i) = vr + vl + npack(i)*(v+mutv)
  end subroutine Update_Rate

  !!$ ----------------------
  !!$ Parallel Using OM
  !!$ ----------------------

  subroutine cell_event_omp(i, kpar)
    use par_zig_mod
    implicit none
    integer, intent(in) :: i, kpar
    integer j, k, m, shift_i
    real u, u1, p0, TGFbeta
    real vr, vl
    type(cell) new_cell
    print *, 'not calling this here'
    read(*,*)
    u = par_uni(kpar)
    u = u*a(i)
    if (npack(i).ne.0) then
       vr = max(0.0, difv*real(npack(i)-npack(i+1)))
       vl = max(0.0, difv*real(npack(i)-npack(i-1)))
    else
       vr = 0.0
       vl = 0.0
    end if
    do j = 1, npack(i)
       u = u - v
       if ( u < 0 ) then
          if ( (cmat(i,j)%type .eq. 1) ) then
             TGFbeta = 0.0
             do k = -brange, brange
                shift_i = k + i
                if ( shift_i .le. 0 ) then
                   shift_i = shift_i + L
                else if ( shift_i > L ) then
                   shift_i = shift_i - L
                end if
                TGFbeta = TGFbeta + TDC(shift_i)*exp(-real(abs(k))/brange)
             end do
             p0 = prelax + (1.0-2.0*prelax) / (1.0 + fdgain1*TGFbeta)
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             u1 = par_uni(kpar)
             if ( u1 < p0 ) then
                ! SC -> 2SC
                cmat(i,j+1) = cmat(i,j)
             else
                ! SC -> 2TAC
                cmat(i,j)%type = 2
                cmat(i,j+1) = cmat(i,j)
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 2 ) then
             ! division
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             u1 = par_uni(kpar)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i, j+1) = cmat(i,j)
             else
                ! TAC -> 2TDC
                cmat(i, j)%type = 3
                cmat(i, j+1) = cmat(i,j)
                TDC(i) = TDC(i) + 2
             end if
             npack(i) = npack(i) + 1
          else if ( cmat(i,j)%type .eq. 3 ) then
             ! death
             do k=j, H-1
                cmat(i, k) = cmat(i, k+1)
             end do
             TDC(i) = TDC(i) - 1
             npack(i) = npack(i) - 1
          else if ( (cmat(i,j)%type .eq. 4) ) then
             p0 = 1.0
             do k=H, j+2, -1
                cmat(i, k) = cmat(i, k-1)
             end do
             cmat(i,j+1) = cmat(i,j)
             npack(i) = npack(i) + 1
          else
             ! do nothing
          end if
          return
       end if
    end do
    u = u - vr
    if ( u < 0 ) then
       u1 = par_uni(kpar)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 4', j, u1
          read(*,*)
       end if

       if ( (cmat(i,j)%type .eq. 1) .or. (cmat(i,j)%type .eq. 4)) then
          u1 = par_uni(kpar)
          if (u1 < scstick) then
             return
          end if
       end if

       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i + 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1
       if (new_cell%type .eq. 3) then
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       end if
       return
    end if
    u = u - vl
    if ( u < 0 ) then
       u1 = par_uni(kpar)
       j = ceiling(u1*npack(i))
       if (j < 1 .or. j>npack(i)) then
          print *, 'error 5', j, u1, npack(i)
          read(*,*)
       end if
       !print *, 'move left at height j', i, j
       if ( (cmat(i,j)%type .eq. 1) .or. (cmat(i,j)%type .eq. 4)) then
          u1 = par_uni(kpar)
          if (u1 < scstick) then
             return
          end if
       end if
       new_cell = cmat(i, j)
       do k=j, H-1
          cmat(i, k) = cmat(i, k+1)
       end do
       npack(i) = npack(i) - 1

       m = i - 1
       if ( cmat(m, j)%type .eq. 0 ) then
          cmat(m, npack(m)+1) = new_cell
       else
          do k=npack(m)+1, j+1, -1
             cmat(m, k) = cmat(m, k-1)
          end do
          cmat(m, j) = new_cell
       end if
       npack(m) = npack(m) + 1
       if (new_cell%type .eq. 3) then
          TDC(i) = TDC(i) - 1
          TDC(m) = TDC(m) + 1
       end if
       return
    else
       print *, 'event at', i
       print *, 'npack', npack(i-1:i+1)
       print *, 'cell', cmat(i, 1:npack(i))%type
       print *, 'u', u
       print *, 'vl', vl
       print *, 'a', a(i)
       write(*,*) 'error 2'
       read(*,*)
    end if
  end subroutine cell_event_omp


  subroutine Next_Reaction_omp(k, tau, ilow, iup)
    implicit none
    integer, intent(out) :: k
    integer, intent(in) :: ilow, iup
    real, intent(out) :: tau
    real tau_temp
    integer i

    tau = huge(0.0)
    k = 0
    do i = ilow, iup
       if ( a(i) > 0.0 ) then
          tau_temp = ( NP(i) - NT(i) ) / a(i)
          if ( tau_temp < tau) then
             tau = tau_temp
             k = i
          end if
       else
!          print *, 'no reaction -------------- ', i, k, a(i)
!          read(*,*)
       end if
    end do
!    print *, 'k', k
!    print *, 'a', a
    if ( k .le. 0 ) then
       write(*,*), 'error'
       read(*,*)
    end if
    if ( tau .le. 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, 'NP-NT', NP - NT
       read(*,*)
    end if
  end subroutine Next_Reaction_omp

  subroutine update_nutri(dt)
    implicit none
    real, intent(in) :: dt
    integer i
    real nutri_flag
    Nutri_old(0:L+1) = Nutri(0:L+1)
    do i = 1, L
       if ( TDC(i) < 0 ) then
          print *, i, TDC(i), SC(i), TAC(i), MC(i), npack(i)
          read(*,*)
          nutri_flag = 0.0
       else
          !nutri_flag = 1.0
          nutri_flag = TDC(i)
       end if
       !nutri_flag = 1.0
       !nutri_flag = TDC(i)
       Nutri(i) = Nutri_old(i) + &
            (Nutri_old(i+1)+Nutri_old(i-1)-2.0*Nutri_old(i))*NutriMobility*dt &
            + NutriGrowthRate*nutri_flag*dt &
            - NutriConsumeRate*(SC(i)+MC(i))*dt &
            - NutriDecayRate*Nutri(i)*dt
       Nutri(i) = max(0.0, Nutri(i))
       !Nutri(i) = min(10.0, Nutri(i))
    end do
    Nutri(L+1) = Nutri(1)
    Nutri(0) = Nutri(L)
    Nutri(1) = 1.0
    Nutri(600) = 1.0
    !Nutri(100) = 15.0
    !Nutri(500) = 15.0
    !print *, Nutri(1:3), SC(1:3)
    !read(*,*)
  end subroutine update_nutri

end module setting
