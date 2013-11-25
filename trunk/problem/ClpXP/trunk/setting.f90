module setting
  real :: k1, k2, r
  real :: jp_base, jp_incre
  integer :: tr_flag, iseed, nsample
  namelist /xdata/ k1, k2, r, jp_base, jp_incre, tr_flag, iseed, nsample

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *), 'Control parameters...'
    write(*, '(a20, f10.2)'), 'k1 = ', k1
    write(*, '(a20, f10.2)'), 'k2 = ', k2
    write(*, '(a20, f10.2)'), 'r = ', r
    write(*, '(a20, f10.2)'), 'jp_base = ', jp_base
    write(*, '(a20, f10.2)'), 'jp_incre = ', jp_incre
    write(*, '(a20, i10)'), 'tr_flag = ', tr_flag
    write(*, '(a20, i10)'), 'iseed = ', iseed
    write(*, '(a20, i10)'), 'nsample = ', nsample

    open(9, file="out/control.csv")
    write(9, '(a20, a10)'), 'PARAMETER,', 'VALUE'
    write(9, '(a20, f10.2)'), 'k1,', k1
    write(9, '(a20, f10.2)'), 'k2,', k2
    write(9, '(a20, f10.2)'), 'r,', r
    write(9, '(a20, f10.2)'), 'jp_base,', jp_base
    write(9, '(a20, f10.2)'), 'jp_incre,', jp_incre
    write(9, '(a20, i10.2)'), 'tr_flag,', tr_flag
    write(9, '(a20, i10.2)'), 'iseed,', iseed
    write(9, '(a20, i10.2)'), 'nsample,', nsample

    close(8)
    close(9)
  end subroutine read_xdata

  subroutine getrate(nATP, a)
    implicit none
    real, intent(in) :: nATP
    real, intent(out) :: a(2)
    
    a = 0.0
    if ( nATP .eq. 0 ) then
       a(1) = k1
       a(2) = 0.0
    else if ( nATP .eq. 1 ) then
       a(1) = k1
       a(2) = 0.0
    else if ( nATP .eq. 2 ) then
       a(1) = k2
       a(2) = 1.0/r
    else if ( nATP .eq. 3 ) then
       a(1) = k2
       a(2) = 1.0/r
    else if ( nATP .eq. 4 ) then
       a(1) = 0
       a(2) = 1.0/r
    else
       print *, 'nATP error'
       read(*,*)
    end if

  end subroutine getrate

end module setting
