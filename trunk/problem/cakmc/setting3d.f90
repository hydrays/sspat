module setting
  integer Lbox
  integer Lbox2
  integer H
  integer :: brange, thisrandseed
  real :: tend, p1, v, difv, mutv, tm
  real :: fdgain1, scstick, prelax, tpinc
  real :: bd10
  integer :: useomp, is64bit
  real :: timestep
  integer :: npar

  namelist /xdata/ Lbox, H, brange, tend, p1, v, difv, mutv, &
       fdgain1, scstick, prelax, thisrandseed, tpinc, tm

  namelist /xdataomp/ useomp, is64bit, timestep, npar

  type cell
     integer type
     real gene1
     real gene2
     real gene3
     real gene4
  end type cell
  type(cell), allocatable :: cmat(:,:,:)
  real, allocatable :: a(:,:), NT(:,:), NP(:,:)
  integer, allocatable :: npack(:, :), TDC(:,:), SC(:,:), TAC(:,:)
  real, allocatable ::  v_diff(:,:,:), geneinfo(:,:,:)
  
contains
  subroutine read_xdata()
    implicit none
    open(8, file="control3d.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)
    read(8, nml=xdataomp)
    
    bd10 = 10.0/real(brange)
    Lbox2 = Lbox*Lbox

    write(*, *), 'Control parameters...'
    write(*, '(a20, i10)'), 'Lbox = ', Lbox
    write(*, '(a20, i10)'), 'H = ', H
    write(*, '(a20, i10)'), 'brange = ', brange
    write(*, '(a20, i10)'), 'thisrandseed = ', thisrandseed
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
    
    if (useomp.eq.1) then
       write(*, '(a)'), 'OpenMP parallel in use!'
       write(*, '(a20, I10)'), 'is64bit = ', is64bit
       write(*, '(a20, f10.2)'), 'timestep = ', timestep
       write(*, '(a20, I10)'), 'npar = ', npar
    end if

    open(9, file="out3d/control3d.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)'), 'Lbox,', Lbox
    write(9, '(a20, i10)'), 'H,', H
    write(9, '(a20, i10)'), 'brange,', brange
    write(9, '(a20, i10)'), 'thisrandseed,', thisrandseed
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
    close(8)
    close(9)
  end subroutine read_xdata

  subroutine init_cell_pool()
    use par_zig_mod
    implicit none
    integer temp_num
    real u
    integer i, j, k

    allocate(cmat(0:Lbox+1,0:Lbox+1,H))
    allocate(a(1:Lbox,1:Lbox))
    allocate(NT(1:Lbox,1:Lbox))
    allocate(NP(1:Lbox,1:Lbox))
    allocate(npack(0:Lbox+1,0:Lbox+1))
    allocate(TDC(0:Lbox+1,0:Lbox+1))
    allocate(SC(0:Lbox+1,0:Lbox+1))
    allocate(TAC(0:Lbox+1,0:Lbox+1))
    allocate(v_diff(1:Lbox, 1:Lbox, 4))
    allocate(geneinfo(0:Lbox+1,0:Lbox+1,4))

    cmat(1:Lbox, 1:Lbox, 1)%type = 1
    cmat(1:Lbox, 1:Lbox, 2:3)%type = 2
    cmat(1:Lbox, 1:Lbox, 3:4)%type = 3
    do i = 1, Lbox
       do j = 1, Lbox
          u = par_uni(0)
          cmat(i,j, 1)%gene1 = scstick!min(u, 0.99)
          cmat(i,j, 1)%gene2 = prelax
          cmat(i,j, 1)%gene3 = fdgain1! + 0.01*(u-0.5)
          cmat(i,j, 1)%gene4 = 10!+2.0*(u-0.5)
       end do
    end do
    cmat(0,1:Lbox,:) = cmat(Lbox,1:Lbox,:)
    cmat(Lbox+1,1:Lbox,:) = cmat(1,1:Lbox,:)
    cmat(1:Lbox,0,:) = cmat(1:Lbox,Lbox,:)
    cmat(1:Lbox,Lbox+1,:) = cmat(1:Lbox,1,:)

    npack = 0
    TDC = 0
    do i = 1, Lbox
       do j = 1, Lbox
          do k = 1, H
             if ( cmat(i,j,k)%type .ne. 0 ) then
                npack(i,j) = npack(i,j) + 1
             end if
             if ( cmat(i,j,k)%type .eq. 3 ) then
                TDC(i,j) = TDC(i,j) + 1
             end if
          end do
       end do
    end do
    npack(0,1:Lbox) = npack(Lbox,1:Lbox)
    npack(Lbox+1,1:Lbox) = npack(1,1:Lbox)
    npack(1:Lbox,0) = npack(1:Lbox,Lbox)
    npack(1:Lbox,Lbox+1) = npack(1:Lbox,1)
    TDC(0,1:Lbox) = TDC(Lbox,1:Lbox)
    TDC(Lbox+1,1:Lbox) = TDC(1,1:Lbox)
    TDC(1:Lbox,0) = TDC(1:Lbox,Lbox)
    TDC(1:Lbox,Lbox+1) = TDC(1:Lbox,1)

    NP = 0.0
    NT = 0.0
    do i = 1, Lbox
       do j = 1, Lbox
          call random_number(u)
          u = -log(u)
          call Update_Rate(i,j)
          NP(i,j) = u
       end do
    end do

  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j, k

    WRITE(filename,'(A9,I5.5,A4)') './out3d/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, Lbox
       do j = 1, Lbox
          geneinfo(i,j,:) = 0.0
          do k = 1, H
!!$             if ( cmat(i,j,k)%type .eq. 1 ) then
!!$                geneinfo(i,j,1) = geneinfo(i,j,1) + cmat(i,j,k)%gene1
!!$                geneinfo(i,j,2) = geneinfo(i,j,2) + cmat(i,j,k)%gene2
!!$                geneinfo(i,j,3) = geneinfo(i,j,3) + cmat(i,j,k)%gene3
!!$                geneinfo(i,j,4) = geneinfo(i,j,4) + cmat(i,j,k)%gene4
!!$             end if
             write(11, '(I5)', advance="no"), cmat(i,j,k)%type
          end do
!!$          if ( SC(i,j).eq.0 ) then
!!$             geneinfo(i,j,:) = -1.0
!!$          else
!!$             geneinfo(i,j,:) = geneinfo(i,j,:) / SC(i,j)
!!$          end if
!!$          write(11, '(f14.7)', advance="no"), geneinfo(i,j,1)
!!$          write(11, '(f14.7)', advance="no"), geneinfo(i,j,3)
          write(11, *)
       end do
    end do
    close(11)
  end subroutine output_to_file

  subroutine cell_event(i, j, kpar)
    use par_zig_mod
    implicit none
    integer, intent(in) :: i, j, kpar
    integer len, k, k2, mi, mj,  shift_i, shift_j
    real u, u1, p0, TGFbeta
    type(cell) new_cell
    u = par_uni(kpar)
    u = u*a(i, j)

!!$    if (kpar .eq. 0) then
!!$       print *, 'event happen at ', kpar, i, j
!!$       print *, 'u', u
!!$       print *, 'a(i,j)', a(i,j)
!!$    print *, 'npack(i)', npack(i-1:i+1)
!!$    print *, 'vr, vl', vr, vl
!!$    end if

    u = u - v_diff(i, j, 1)
    if ( u < 0 ) then
       mi = i+1
       mj = j
    else
       u = u - v_diff(i, j, 2)
       if ( u < 0 ) then
          mi = i-1
          mj = j
       else
          u = u - v_diff(i, j, 3)
          if ( u < 0 ) then
             mi = i
             mj = j+1
          else
             u = u - v_diff(i, j, 4)
             if ( u < 0 ) then
                mi = i
                mj = j-1
             end if
          end if
       end if
    end if
    if ( u < 0 ) then
       u1 = par_uni(kpar)
       len = ceiling(u1*npack(i, j))
       if (len < 1 .or. len >npack(i,j)) then
          print *, 'error 4', len, u1
          read(*,*)
       end if

       if ( cmat(i,j,len)%type .eq. 1 ) then
          u1 = par_uni(kpar)
          if ( u1 < cmat(i,j,len)%gene1 ) then
             return
          end if
       end if

       new_cell = cmat(i, j, len)
       do k=len, H-1
          cmat(i,j,k) = cmat(i,j,k+1)
       end do
       npack(i,j) = npack(i,j) - 1

       if ( cmat(mi,mj,len)%type .eq. 0 ) then
          cmat(mi,mj,npack(mi,mj)+1) = new_cell
       else
          do k=npack(mi,mj)+1, len+1, -1
             cmat(mi,mj,k) = cmat(mi,mj,k-1)
          end do
          cmat(mi,mj,len) = new_cell
       end if
       npack(mi,mj) = npack(mi,mj) + 1
       if (new_cell%type .eq. 3) then
          TDC(i,j) = TDC(i,j) - 1
          TDC(mi,mj) = TDC(mi,mj) + 1
       end if
       return
    end if

    do len = 1, npack(i,j)
       u = u - v
       if ( u < 0 ) then
          if ( cmat(i,j,len)%type .eq. 1 ) then
             TGFbeta = 0.0
             do k = -brange, brange
                do k2 = -brange, brange
                   shift_i = k + i
                   shift_j = k2 + j
                   if ( shift_i .le. 0 ) then
                      shift_i = shift_i + Lbox
                   else if ( shift_i > Lbox ) then
                      shift_i = shift_i - Lbox
                   end if
                   if ( shift_j .le. 0 ) then
                      shift_j = shift_j + Lbox
                   else if ( shift_j > Lbox ) then
                      shift_j = shift_j - Lbox
                   end if
                   TGFbeta = TGFbeta + bd10*TDC(shift_i,shift_j)*&
                        exp(- sqrt( (real(abs(k))/brange)**2 + &
                        (real(abs(k2))/brange)**2 ) )
                end do
             end do
             p0 = cmat(i,j,len)%gene2 + (1.0 - 2.0*cmat(i,j,len)%gene2) &
                  / (1.0 + cmat(i,j,len)%gene3*TGFbeta)
             !p0 = 0.2 + 0.6 / (1.0 + 0.005*TGFbeta)

             !print *, 'p0', p0
             ! division
             do k=H, len+2, -1
                cmat(i,j,k) = cmat(i,j,k-1)
             end do
             u1 = par_uni(kpar)

             !write(*, *), kpar, u1

             if ( u1 < p0 ) then
                ! SC -> 2SC
                cmat(i,j,len+1) = cmat(i,j,len)
             else
                ! SC -> 2TAC
                cmat(i,j,len)%type = 2
                cmat(i,j,len+1) = cmat(i,j,len) 
             end if
             npack(i,j) = npack(i,j) + 1
          else if ( cmat(i,j,len)%type .eq. 2 ) then
             ! division
             do k=H, len+2, -1
                cmat(i,j,k) = cmat(i,j,k-1)
             end do
             u1 = par_uni(kpar)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i,j,len+1) = cmat(i,j,len)
             else
                ! TAC -> 2TDC
                cmat(i,j,len)%type = 3
                cmat(i,j,len+1) = cmat(i,j,len)
                TDC(i,j) = TDC(i,j) + 2
             end if
             npack(i,j) = npack(i,j) + 1
          else if ( cmat(i,j,len)%type .eq. 3 ) then
             ! death
             do k=len, H-1
                cmat(i,j,k) = cmat(i,j,k+1)
             end do
             TDC(i,j) = TDC(i,j) - 1
             npack(i,j) = npack(i,j) - 1
          else
             ! do nothing
          end if
          return
       end if
    end do
    if ( mutv .ne. 0.0 ) then
       do len = 1, npack(i,j)
          u = u - mutv
          if ( u < 0 ) then
             if ( cmat(i,j,len)%type .eq. 1 ) then
                ! mutation
                u1 = par_uni(kpar)
                cmat(i,j,len)%gene1 = cmat(i,j,len)%gene1 + (u1-0.5)*0.1
                if ( cmat(i,j,len)%gene1 > 0.995) then
                   cmat(i,j,len)%gene1 = 0.995
                else if ( cmat(i,j,len)%gene1 < 0.0) then
                   cmat(i,j,len)%gene1 = 0.0
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
    print *, "u", u
    print *, "a", a(i,j)
    print *, 'i,j', i, j
    stop
  end subroutine cell_event

  subroutine cell_stat(t)
    implicit none

    real, intent(in) :: t
    integer i, j, k

    SC = 0
    TAC = 0
    do i = 1, Lbox
       do j = 1, Lbox
          do k = 1, H
             if (cmat(i,j,k)%type.eq.1) then
                SC(i,j) = SC(i,j) + 1
             elseif (cmat(i,j,k)%type.eq.2) then
                TAC(i,j) = TAC(i,j) + 1
             end if
          end do
       end do
    end do

    write(*, '(5(F10.2))'), t, real(sum(SC))/Lbox2, &
         real(sum(tac))/Lbox2, real(sum(tdc))/Lbox2
    write(100, '(5(F10.2))'), t, real(sum(SC))/Lbox2, &
         real(sum(tac))/Lbox2, real(sum(tdc))/Lbox2
  end subroutine cell_stat

  subroutine Perodic_BC(i, j)
    implicit none
    integer, intent(in) :: i, j

    if ( i.eq.1 ) then
       if ( j.eq.1 ) then !upper-left coner
          cmat(Lbox:Lbox+1, 1:2, :) = cmat(0:1, 1:2, :) 
          npack(Lbox:Lbox+1, 1:2) = npack(0:1,1:2)
          TDC(Lbox:Lbox+1,1:2) = TDC(0:1,1:2)

          cmat(1:2, Lbox:Lbox+1, :) = cmat(1:2, 0:1, :) 
          npack(1:2, Lbox:Lbox+1) = npack(1:2, 0:1)
          TDC(1:2, Lbox:Lbox+1) = TDC(1:2, 0:1)
       else if ( j.eq.Lbox ) then !down-left coner
          cmat(Lbox:Lbox+1, Lbox-1:Lbox, :) = cmat(0:1, Lbox-1:Lbox, :) 
          npack(Lbox:Lbox+1, Lbox-1:Lbox) = npack(0:1,Lbox-1:Lbox)
          TDC(Lbox:Lbox+1,Lbox-1:Lbox) = TDC(0:1,Lbox-1:Lbox)

          cmat(1:2, 0:1, :) = cmat(1:2, Lbox:Lbox+1, :) 
          npack(1:2, 0:1) = npack(1:2, Lbox:Lbox+1)
          TDC(1:2, 0:1) = TDC(1:2, Lbox:Lbox+1)
       else ! left edge
          cmat(Lbox, j-1:j+1, :) = cmat(0, j-1:j+1, :) 
          cmat(Lbox+1, j-1:j+1, :) = cmat(1, j-1:j+1, :)
          npack(Lbox, j-1:j+1) = npack(0,j-1:j+1)
          npack(Lbox+1,j-1:j+1) = npack(1,j-1:j+1)
          TDC(Lbox,j-1:j+1) = TDC(0,j-1:j+1)
          TDC(Lbox+1,j-1:j+1) = TDC(1,j-1:j+1)
       end if
    end if
    if ( i .eq. 2 ) then
       cmat(Lbox+1,j-1:j+1, :) = cmat(1, j-1:j+1,:)
       npack(Lbox+1,j-1:j+1) = npack(1,j-1:j+1)
       TDC(Lbox+1,j-1:j+1) = TDC(1,j-1:j+1)
    end if

    if ( i .eq. Lbox ) then
       if ( j.eq.1 ) then !upper-right coner
          cmat(0:1, 1:2, :) = cmat(Lbox:Lbox+1, 1:2, :) 
          npack(0:1, 1:2) = npack(Lbox:Lbox+1,1:2)
          TDC(0:1, 1:2) = TDC(Lbox:Lbox+1,1:2)

          cmat(Lbox-1:Lbox, Lbox:Lbox+1, :) = cmat(Lbox-1:Lbox, 0:1, :) 
          npack(Lbox-1:Lbox, Lbox:Lbox+1) = npack(Lbox-1:Lbox, 0:1)
          TDC(Lbox-1:Lbox, Lbox:Lbox+1) = TDC(Lbox-1:Lbox, 0:1)
       else if ( j.eq.Lbox ) then !down-right coner
          cmat(0:1, Lbox-1:Lbox, :) = cmat(Lbox:Lbox+1, Lbox-1:Lbox, :) 
          npack(0:1, Lbox-1:Lbox) = npack(Lbox:Lbox+1,Lbox-1:Lbox)
          TDC(0:1,Lbox-1:Lbox) = TDC(Lbox:Lbox+1,Lbox-1:Lbox)

          cmat(Lbox-1:Lbox, 0:1, :) = cmat(Lbox-1:Lbox, Lbox:Lbox+1, :) 
          npack(Lbox-1:Lbox, 0:1) = npack(Lbox-1:Lbox, Lbox:Lbox+1)
          TDC(Lbox-1:Lbox, 0:1) = TDC(Lbox-1:Lbox, Lbox:Lbox+1)
       else ! right edge
          cmat(0:1, j-1:j+1,:) = cmat(Lbox:Lbox+1, j-1:j+1,:)
          npack(0:1,j-1:j+1) = npack(Lbox:Lbox+1,j-1:j+1)
          TDC(0:1,j-1:j+1) = TDC(Lbox:Lbox+1,j-1:j+1)
       end if
    end if
    if ( i .eq. Lbox-1 ) then
       cmat(0, j-1:j+1,:) = cmat(Lbox, j-1:j+1,:)
       npack(0,j-1:j+1) = npack(Lbox,j-1:j+1)
       TDC(0,j-1:j+1) = TDC(Lbox,j-1:j+1)
    end if


    if ( j .eq. 1 ) then
       if ( i .eq. 1 ) then
          ! do nothing
       else if ( i.eq. Lbox) then
          ! do nothing
       else
          cmat(i-1:i+1,Lbox:Lbox+1, :) = cmat(i-1:i+1,0:1, :) 
          npack(i-1:i+1,Lbox:Lbox+1) = npack(i-1:i+1,0:1)
          TDC(i-1:i+1,Lbox:Lbox+1) = TDC(i-1:i+1,0:1)
       end if
    end if
    if ( j .eq. 2 ) then
       cmat(i-1:i+1,Lbox+1, :) = cmat(i-1:i+1,1, :)
       npack(i-1:i+1,Lbox+1) = npack(i-1:i+1,1)
       TDC(i-1:i+1,Lbox+1) = TDC(i-1:i+1,1)
    end if

    if ( j .eq. Lbox ) then
       if ( i .eq. 1 ) then
          ! do nothing
       else if ( i.eq. Lbox) then
          ! do nothing
       else
          cmat(i-1:i+1,0:1, :) = cmat(i-1:i+1,Lbox:Lbox+1, :)
          npack(i-1:i+1,0:1) = npack(i-1:i+1,Lbox:Lbox+1)
          TDC(i-1:i+1,0:1) = TDC(i-1:i+1,Lbox:Lbox+1)
       end if
    end if
    if ( j .eq. Lbox-1 ) then
       cmat(i-1:i+1,0, :) = cmat(i-1:i+1,Lbox, :)
       npack(i-1:i+1,0) = npack(i-1:i+1,Lbox)
       TDC(i-1:i+1,0) = TDC(i-1:i+1,Lbox)
    end if
  end subroutine Perodic_BC

  subroutine Next_Reaction(k1, k2, tau, ilow, iup, jlow, jup)
    implicit none
    integer, intent(out) :: k1, k2
    integer, intent(in) :: ilow, iup, jlow, jup
    real, intent(out) :: tau
    real tau_temp
    integer i, j

    tau = huge(0.0)
    k1 = 0
    k2 = 0
    do j = jlow, jup
       do i = ilow, iup
          if ( a(i,j) > 0.0 ) then
             tau_temp = ( NP(i,j) - NT(i,j) ) / a(i,j)
             if ( tau_temp < tau) then
                tau = tau_temp
                k1 = i
                k2 = j
             end if
          else
          !          print *, 'no reaction -------------- ', i, k, a(i)
          !          read(*,*)
          end if
       end do
    end do
    !    print *, 'k', k
    !    print *, 'a', a
    if ( k1 < 0 .or. k2 < 0 ) then
       write(*,*), 'error'
       read(*,*)
    end if
    if ( tau < 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, 'NP-NT', NP - NT
       read(*,*)
    end if
  end subroutine Next_Reaction

  subroutine Update_Rate(i,j)
    implicit none
    integer, intent(in) :: i,j
    v_diff(i, j, 1) = max(0.0, difv*real(npack(i,j)-npack(i+1,j)))
    v_diff(i, j, 2) = max(0.0, difv*real(npack(i,j)-npack(i-1,j)))
    v_diff(i, j, 3) = max(0.0, difv*real(npack(i,j)-npack(i,j+1)))
    v_diff(i, j, 4) = max(0.0, difv*real(npack(i,j)-npack(i,j-1)))
    a(i,j) = sum(v_diff(i, j, :)) + npack(i,j)*(v+mutv)
  end subroutine Update_Rate

end module setting
