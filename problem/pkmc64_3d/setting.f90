module setting
  integer, parameter :: Lbox = 128
  integer, parameter :: Lbox2 = Lbox*Lbox
  integer, parameter :: H = 100
  real, parameter :: b = 8.0
  real, parameter :: bd10 = 10.0/b
  real, parameter :: tend = 500.0
  real, parameter :: p1 = 0.3
  real, parameter :: v = 1.0
  real, parameter :: D = 50.0
  real, parameter :: mv = 0.0
  type cell
     integer type
     real gene1
     real gene2
     real gene3
     real gene4
  end type cell
  type(cell), allocatable :: cmat(:,:,:)
!  type(cell) cmat(0:Lbox+1,0:Lbox+1,H)
  real, allocatable :: a(:,:), NT(:,:), NP(:,:)
  integer, allocatable :: npack(:, :), TDC(:,:), SC(:,:), TAC(:,:)
  real, allocatable ::  v_diff(:,:,:), geneinfo(:,:,:)
  
contains
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
          !u = par_uni(0)
          cmat(i,j, 1)%gene1 = 0.96!min(u, 0.96)
          cmat(i,j, 1)%gene2 = 0.2
          cmat(i,j, 1)%gene3 = 0.001! + 0.01*(u-0.5)
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

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, Lbox+1
       do j = 1, Lbox+1
          geneinfo(i,j,:) = 0.0
          do k = 1, H
             if ( cmat(i,j,k)%type .eq. 1 ) then
                geneinfo(i,j,1) = geneinfo(i,j,1) + cmat(i,j,k)%gene1
                geneinfo(i,j,2) = geneinfo(i,j,2) + cmat(i,j,k)%gene2
                geneinfo(i,j,3) = geneinfo(i,j,3) + cmat(i,j,k)%gene3
                geneinfo(i,j,4) = geneinfo(i,j,4) + cmat(i,j,k)%gene4
             end if
             write(11, '(I5)', advance="no"), cmat(i,j,k)%type
          end do
          if ( SC(i,j).eq.0 ) then
             geneinfo(i,j,:) = -1.0
          else
             geneinfo(i,j,:) = geneinfo(i,j,:) / SC(i,j)
          end if
          write(11, '(I6)', advance="no"), TDC(i,j)
          write(11, '(I6)', advance="no"), TAC(i,j)
          write(11, '(I6)', advance="no"), SC(i,j)
          write(11, '(I6)', advance="no"), npack(i,j)
          write(11, '(f14.7)', advance="no"), geneinfo(i,j,1)
          write(11, *)
       end do
    end do
    close(11)
  end subroutine output_to_file

  subroutine cell_event(i, j, kpar)
    use par_zig_mod
    implicit none
    integer, intent(in) :: i, j, kpar
    integer l, k, k2, mi, mj,  shift_i, shift_j
    real u, u1, p0, TGFbeta
    type(cell) new_cell
    u = par_uni(kpar)
    u = u*a(i, j)
!!$    print *, 'event happen at ', i
!!$    print *, 'u', u
!!$    print *, 'a(i)', a(i)
!!$    print *, 'npack(i)', npack(i-1:i+1)
!!$    print *, 'vr, vl', vr, vl

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
       l = ceiling(u1*npack(i, j))
       if (l < 1 .or. l>npack(i,j)) then
          print *, 'error 4', l, u1
          read(*,*)
       end if

       if ( cmat(i,j,l)%type .eq. 1 ) then
          u1 = par_uni(kpar)
          if ( u1 < cmat(i,j,l)%gene1 ) then
             return
          end if
       end if

       new_cell = cmat(i, j, l)
       do k=l, H-1
          cmat(i,j,k) = cmat(i,j,k+1)
       end do
       npack(i,j) = npack(i,j) - 1

       if ( cmat(mi,mj,l)%type .eq. 0 ) then
          cmat(mi,mj,npack(mi,mj)+1) = new_cell
       else
          do k=npack(mi,mj)+1, l+1, -1
             cmat(mi,mj,k) = cmat(mi,mj,k-1)
          end do
          cmat(mi,mj,l) = new_cell
       end if
       npack(mi,mj) = npack(mi,mj) + 1
       if (new_cell%type .eq. 3) then
          TDC(i,j) = TDC(i,j) - 1
          TDC(mi,mj) = TDC(mi,mj) + 1
       end if
       return
    end if

    do l = 1, npack(i,j)
       u = u - v
       if ( u < 0 ) then
          if ( cmat(i,j,l)%type .eq. 1 ) then
             TGFbeta = 0.0
             do k = -b, b
                do k2 = -b, b
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
                        exp(- sqrt( (real(abs(k))/b)**2 + (real(abs(k2))/b)**2 ) )
                end do
             end do
             p0 = cmat(i,j,l)%gene2 + (1.0 - 2.0*cmat(i,j,l)%gene2) &
                  / (1.0 + cmat(i,j,l)%gene3*TGFbeta)
             !p0 = 0.2 + 0.6 / (1.0 + 0.005*TGFbeta)

             !print *, 'p0', p0
             ! division
             do k=H, l+2, -1
                cmat(i,j,k) = cmat(i,j,k-1)
             end do
             u1 = par_uni(kpar)
             if ( u1 < p0 ) then
                ! SC -> 2SC
                cmat(i,j,l+1) = cmat(i,j,l)
             else
                ! SC -> 2TAC
                cmat(i,j,l)%type = 2
                cmat(i,j,l+1) = cmat(i,j,l) 
             end if
             npack(i,j) = npack(i,j) + 1
          else if ( cmat(i,j,l)%type .eq. 2 ) then
             ! division
             do k=H, l+2, -1
                cmat(i,j,k) = cmat(i,j,k-1)
             end do
             u1 = par_uni(kpar)
             if ( u1 < p1 ) then
                ! TAC -> 2TAC
                cmat(i,j,l+1) = cmat(i,j,l)
             else
                ! TAC -> 2TDC
                cmat(i,j,l)%type = 3
                cmat(i,j,l+1) = cmat(i,j,l)
                TDC(i,j) = TDC(i,j) + 2
             end if
             npack(i,j) = npack(i,j) + 1
          else if ( cmat(i,j,l)%type .eq. 3 ) then
             ! death
             do k=l, H-1
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
    if ( mv .ne. 0.0 ) then
       do l = 1, npack(i,j)
          u = u - mv
          if ( u < 0 ) then
             if ( cmat(i,j,l)%type .eq. 1 ) then
                ! mutation
                u1 = par_uni(kpar)
                cmat(i,j,l)%gene1 = cmat(i,j,l)%gene1 + (u1-0.5)*0.1
                if ( cmat(i,j,l)%gene1 > 0.995) then
                   cmat(i,j,l)%gene1 = 0.995
                else if ( cmat(i,j,l)%gene1 < 0.0) then
                   cmat(i,j,l)%gene1 = 0.0
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
    if ( k1 .le. 0 .or. k2 .le. 0 ) then
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
    v_diff(i, j, 1) = max(0.0, D*real(npack(i,j)-npack(i+1,j)))
    v_diff(i, j, 2) = max(0.0, D*real(npack(i,j)-npack(i-1,j)))
    v_diff(i, j, 3) = max(0.0, D*real(npack(i,j)-npack(i,j+1)))
    v_diff(i, j, 4) = max(0.0, D*real(npack(i,j)-npack(i,j-1)))
    a(i,j) = sum(v_diff(i, j, :)) + npack(i,j)*(v+mv)
  end subroutine Update_Rate

end module setting
