module setting
  integer Lbox
  integer :: iseed
  real :: tend, dt
  real :: alpha, beta, alpha_max
  integer :: R1, nc
  real :: diff, lambda, gamma
  real :: k_lambda
  integer :: model_type, read_lambda_from_file

  namelist /xdata/ Lbox, tend, dt, alpha, beta, alpha_max, diff, k_lambda, &
       iseed, tpinc, R1, nc, lambda, gamma, model_type, &
       read_lambda_from_file

  type cell
     integer type
     ! type == 0 : empty splot
     ! type == 1 : blood_vessel
     ! type == 2 : liver cell
     ! type == 5 : source (artery)
     ! type == 6 : sink (vein)       
     integer n
     ! number of blood cells in the grid
     real z
     ! size of the vessel
     real phi
     ! direction of the vessel
  end type cell

  type(cell), allocatable :: cmat(:,:)
  real, allocatable :: p(:,:)
  real, allocatable :: a(:,:)

contains
  subroutine read_xdata()
    implicit none
    open(8, file="control.txt", status='OLD', recl=80, delim='APOSTROPHE')
    read(8, nml=xdata)

    write(*, *) 'Control parameters...'
    write(*, '(a20, i10)') 'Lbox = ', Lbox
    write(*, '(a20, i10)') 'iseed = ', iseed
    write(*, '(a20, f10.2)') 'tend = ', tend
    write(*, '(a20, f10.2)') 'dt = ', dt
    write(*, '(a20, f10.2)') 'alpha = ', alpha
    write(*, '(a20, f10.2)') 'beta = ', beta
    write(*, '(a20, f10.2)') 'alpha_max = ', alpha_max
    write(*, '(a20, f10.2)') 'diff = ', diff
    write(*, '(a20, f10.2)') 'tpinc = ', tpinc
    write(*, '(a20, i10)') 'R1 = ', R1
    write(*, '(a20, f10.2)') 'lambda = ', lambda
    write(*, '(a20, f10.2)') 'gamma = ', gamma
    write(*, '(a20, f10.2)') 'k_lambda = ', k_lambda
    write(*, '(a20, i10)') 'nc = ', nc
    write(*, '(a20, i10)') 'model_type = ', model_type
    write(*, '(a20, i10)') 'read_lambda_from_file = ', read_lambda_from_file

    open(9, file="out/control.csv")
    write(9, '(a20, a10)') 'PARAMETER,', 'VALUE'
    write(9, '(a20, i10)') 'Lbox,', Lbox
    write(9, '(a20, i10)') 'iseed,', iseed
    write(9, '(a20, f10.2)') 'tend,', tend
    write(9, '(a20, f10.2)') 'dt,', dt
    write(9, '(a20, f10.2)') 'alpha,', alpha
    write(9, '(a20, f10.2)') 'beta,', beta
    write(9, '(a20, f10.2)') 'alpha_max,', alpha_max
    write(9, '(a20, f10.2)') 'diff,', diff
    write(9, '(a20, f10.2)') 'tpinc,', tpinc
    write(9, '(a20, i10)') 'R1,', R1
    write(9, '(a20, f10.2)') 'lambda,', lambda
    write(9, '(a20, f10.2)') 'gamma,', gamma
    write(9, '(a20, f10.2)') 'k_lambda,', k_lambda
    write(9, '(a20, i10)') 'nc,', nc
    write(9, '(a20, i10)') 'model_type,', model_type
    write(9, '(a20, i10)') 'read_lambda_from_file,', read_lambda_from_file

    close(8)
    close(9)
  end subroutine read_xdata

  subroutine init_cell_pool()
    implicit none
    real x, y
    real u1, u2
    integer i, j, ishift, jshift
    integer curb

    write(*, '(A)', advance='no') 'Initialize...'
    allocate(cmat(1:Lbox,1:Lbox))
    allocate(p(1:Lbox,1:Lbox))
    allocate(a(1:Lbox,1:Lbox))

    do i = 1, Lbox
      y = i - Lbox/2      
      do j = 1, Lbox
         x = j - Lbox/2
         cmat(i,j)%type = 1
         cmat(i,j)%phi = atan2(x, y)
         cmat(i,j)%z = abs(x-y)
         cmat(i,j)%n = 0
      end do
    end do

    cmat(Lbox/2, Lbox/2)%type = 5
    cmat(Lbox-10, Lbox-10)%type = 6
    cmat(1+10, 1+10)%type = 6
    cmat(1+10, Lbox-10)%type = 6
    cmat(Lbox-10, 1+10)%type = 6
    cmat(1, 1:Lbox)%type = 6
    cmat(Lbox, 1:Lbox)%type = 6
    cmat(1:Lbox, 1)%type = 6
    cmat(1:Lbox, Lbox)%type = 6

    call update_n()
    write(*, *) 'Done.'
  end subroutine init_cell_pool

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i, j
    real x, y

    WRITE(filename,'(A7,I5.5,A4)') './out/a', index, '.dat'
    open (unit = 11, file=filename, action="write")
    write(11, '(A4, A2, A4, A2, A4, A2, A4, A2, A4, A2, A4)') &
    'x', ', ', 'y', ', ', 'type', ', ', &
    'z', ', ', 'phi', ', ', 'n'

    do i = 1, Lbox
       y = i - Lbox/2      
       do j = 1, Lbox
          x = j - Lbox/2
          write(11, '(F12.4, A2, F12.4, A2, I5, A2, F12.4, A2, F12.4, A2, I12)') &
            x, ', ', y, ', ', cmat(i,j)%type, ', ', &
            cmat(i,j)%z, ', ', cmat(i,j)%phi, ', ', cmat(i,j)%n
       end do
    end do
    close(11)
  end subroutine output_to_file

  subroutine cell_event(i, j)
    implicit none
    integer, intent(in) :: i, j
  end subroutine cell_event

  subroutine cell_stat()
    implicit none
  end subroutine cell_stat

  subroutine update_n()
    implicit none
    integer i, j
    real p(9)
    integer ntot
    integer ix(9)
    real ptot, prob
    integer icat
    do i = 1, Lbox
       do j = 1, Lbox
          if (cmat(i, j)%type == 5) then
             cmat(i, j)%n = cmat(i, j)%n + 100
          else if (cmat(i, j)%type == 6) then
             cmat(i, j)%n = cmat(i, j)%n - min(100, cmat(i, j)%n)
          endif
       end do
    end do
    
    do i = 2, Lbox-1
      do j = 2, Lbox-1
        if (cmat(i, j)%n > 0) then
            call get_p(p)
            !
            !  Initialize variables.
            !
            ntot = cmat(i, j)%n
            ptot = 1.0
            do icat = 1, 9
               ix(icat) = 0
            end do
            !
            !  Generate the observation.
            !
            do icat = 1, 8
               prob = p(icat) / ptot
               ix(icat) = ignbin ( ntot, prob )
               ntot = ntot - ix(icat)
               if ( ntot <= 0 ) then
                  exit
               end if
               ptot = ptot - p(icat)
            end do
            if (ntot > 0) then
               ix(9) = ntot
            endif
            ! check multinomial distribution
            if ( sum(ix) .ne. cmat(i, j)%n ) then
               print *, 'multinomial wrong'
               exit
            endif
            if ( ntot < 0 ) then
               print *, 'multinomial wrong'
               exit
            endif
            cmat(i, j)%n = ix(1)
            cmat(i, j+1)%n = cmat(i, j+1)%n + ix(2)
            cmat(i+1, j+1)%n = cmat(i+1, j+1)%n + ix(3)
            cmat(i+1, j)%n = cmat(i+1, j)%n + ix(4)
            cmat(i+1, j-1)%n = cmat(i+1, j-1)%n + ix(5)
            cmat(i, j-1)%n = cmat(i, j-1)%n + ix(6)
            cmat(i-1, j-1)%n = cmat(i-1, j-1)%n + ix(7)
            cmat(i-1, j)%n = cmat(i-1, j)%n + ix(8)
            cmat(i-1, j+1)%n = cmat(i-1, j+1)%n + ix(9)
         end if
      end do
   end do    
  end subroutine update_n

  subroutine get_p(i, j, p)
    implicit none
    integer, intent(in) :: i, j
    real, intent(out) :: p(9)
    real s0 = 1.0
    p(1) = s0 + cmat(i, j)%z * cmat(i, j)%z
    p(2) = s0 + cmat(i, j)%z * cmat(i, j+1)%z   * cos(cmat(i, j)%phi - cmat(i, j+1)%phi)
    p(3) = s0 + cmat(i, j)%z * cmat(i+1, j+1)%z * cos(cmat(i, j)%phi - cmat(i+1, j+1)%phi)
    p(4) = s0 + cmat(i, j)%z * cmat(i+1, j)%z   * cos(cmat(i, j)%phi - cmat(i+1, j)%phi)
    p(5) = s0 + cmat(i, j)%z * cmat(i+1, j-1)%z * cos(cmat(i, j)%phi - cmat(i+1, j-1)%phi)
    p(6) = s0 + cmat(i, j)%z * cmat(i, j-1)%z   * cos(cmat(i, j)%phi - cmat(i, j-1)%phi)
    p(7) = s0 + cmat(i, j)%z * cmat(i-1, j-1)%z * cos(cmat(i, j)%phi - cmat(i-1, j-1)%phi)
    p(8) = s0 + cmat(i, j)%z * cmat(i-1, j)%z   * cos(cmat(i, j)%phi - cmat(i-1, j)%phi)
    p(9) = s0 + cmat(i, j)%z * cmat(i-1, j+1)%z * cos(cmat(i, j)%phi - cmat(i-1, j+1)%phi)

    cmat(i, j+1)%n = cmat(i, j+1)%n + ix(2)
    cmat(i+1, j+1)%n = cmat(i+1, j+1)%n + ix(3)
    cmat(i+1, j)%n = cmat(i+1, j)%n + ix(4)
    cmat(i+1, j-1)%n = cmat(i+1, j-1)%n + ix(5)
    cmat(i, j-1)%n = cmat(i, j-1)%n + ix(6)
    cmat(i-1, j-1)%n = cmat(i-1, j-1)%n + ix(7)
    cmat(i-1, j)%n = cmat(i-1, j)%n + ix(8)
    cmat(i-1, j+1)%n = cmat(i-1, j+1)%n + ix(9)
    if ( abs(sum(p) - 1) > 1e-12 ) then
       stop 'wrong probabilites'
    end if
  end subroutine get_p
  
  function ignbin ( n, pp )

   !*****************************************************************************80
   !
   !! IGNBIN generates a binomial random deviate.
   !
   !  Discussion:
   !
   !    This procedure generates a single random deviate from a binomial
   !    distribution whose number of trials is N and whose
   !    probability of an event in each trial is P.
   !
   !    The previous version of this program relied on the assumption that
   !    local memory would be preserved between calls.  It set up data
   !    one time to be preserved for use over multiple calls.  In the
   !    interests of portability, this assumption has been removed, and
   !    the "setup" data is recomputed on every call.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    31 March 2013
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Barry Brown, James Lovato.
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Voratas Kachitvichyanukul, Bruce Schmeiser,
   !    Binomial Random Variate Generation,
   !    Communications of the ACM,
   !    Volume 31, Number 2, February 1988, pages 216-222.
   !
   !  Parameters:
   !
   !    Input, integer N, the number of binomial trials, from which a
   !    random deviate will be generated.
   !    0 < N.
   !
   !    Input, real PP, the probability of an event in each trial of
   !    the binomial distribution from which a random deviate is to be generated.
   !    0.0 < PP < 1.0.
   !
   !    Output, integer IGNBIN, a random deviate from the
   !    distribution.
   !
     implicit none
   
     real al
     real alv
     real amaxp
     real c
     real f
     real f1
     real f2
     real ffm
     real fm
     real g
     integer i
     integer ignbin
     integer ix
     integer ix1
     integer k
     integer m
     integer mp
     real pp
     integer n
     real p
     real p1
     real p2
     real p3
     real p4
     real q
     real qn
     real r
     real r4_uni_01
     real t
     real u
     real v
     real w
     real w2
     real x
     real x1
     real x2
     real xl
     real xll
     real xlr
     real xm
     real xnp
     real xnpq
     real xr
     real ynorm
     real z
     real z2
   
     if ( pp <= 0.0E+00 .or. 1.0E+00 <= pp ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'IGNBIN - Fatal error!'
       write ( *, '(a)' ) '  PP is out of range.'
       stop 1
     end if
   
     p = min ( pp, 1.0E+00 - pp )
     q = 1.0E+00 - p
     xnp = real ( n, kind = 4 ) * p
   
     if ( xnp < 30.0E+00 ) then
   
       qn = q ** n
       r = p / q
       g = r * real ( n + 1, kind = 4 )
   
       do
   
         ix = 0
         f = qn
         call random_number(u)
   
         do
   
           if ( u < f ) then
             if ( 0.5E+00 < pp ) then
               ix = n - ix
             end if
             ignbin = ix
             return
           end if
   
           if ( 110 < ix ) then
             exit
           end if
   
           u = u - f
           ix = ix + 1
           f = f * ( g / real ( ix, kind = 4 ) - r )
   
         end do
   
       end do
   
     end if
   
     ffm = xnp + p
     m = int ( ffm )
     fm = m
     xnpq = xnp * q
     p1 = int ( 2.195E+00 * sqrt ( xnpq ) - 4.6E+00 * q ) + 0.5E+00
     xm = fm + 0.5E+00
     xl = xm - p1
     xr = xm + p1
     c = 0.134E+00 + 20.5E+00 / ( 15.3E+00 + fm )
     al = ( ffm - xl ) / ( ffm - xl * p )
     xll = al * ( 1.0E+00 + 0.5E+00 * al )
     al = ( xr - ffm ) / ( xr * q )
     xlr = al * ( 1.0E+00 + 0.5E+00 * al )
     p2 = p1 * ( 1.0E+00 + c + c )
     p3 = p2 + c / xll
     p4 = p3 + c / xlr
   !
   !  Generate a variate.
   !
     do
   
       call random_number(u)
       u = u * p4
       call random_number(v)
   !
   !  Triangle
   !
       if ( u < p1 ) then
         ix = int ( xm - p1 * v + u )
         if ( 0.5E+00 < pp ) then
           ix = n - ix
         end if
         ignbin = ix
         return
       end if
   !
   !  Parallelogram
   !
       if ( u <= p2 ) then
   
         x = xl + ( u - p1 ) / c
         v = v * c + 1.0E+00 - abs ( xm - x ) / p1
   
         if ( v <= 0.0E+00 .or. 1.0E+00 < v ) then
           cycle
         end if
   
         ix = int ( x )
   
       else if ( u <= p3 ) then
   
         ix = int ( xl + log ( v ) / xll )
         if ( ix < 0 ) then
           cycle
         end if
         v = v * ( u - p2 ) * xll
   
       else
   
         ix = int ( xr - log ( v ) / xlr )
         if ( n < ix ) then
           cycle
         end if
         v = v * ( u - p3 ) * xlr
   
       end if
   
       k = abs ( ix - m )
   
       if ( k <= 20 .or. xnpq / 2.0 - 1.0 <= k ) then
   
         f = 1.0E+00
         r = p / q
         g = ( n + 1 ) * r
   
         if ( m < ix ) then
           mp = m + 1
           do i = m + 1, ix
             f = f * ( g / i - r )
           end do
         else if ( ix < m ) then
           ix1 = ix + 1
           do i = ix + 1, m
             f = f / ( g / real ( i, kind = 4 ) - r )
           end do
         end if
   
         if ( v <= f ) then
           if ( 0.5E+00 < pp ) then
             ix = n - ix
           end if
           ignbin = ix
           return
         end if
   
       else
   
         amaxp = ( k / xnpq ) * ( ( k * ( k / 3.0E+00 &
           + 0.625E+00 ) + 0.1666666666666E+00 ) / xnpq + 0.5E+00 )
         ynorm = - real ( k * k, kind = 4 ) / ( 2.0E+00 * xnpq )
         alv = log ( v )
   
         if ( alv < ynorm - amaxp ) then
           if ( 0.5E+00 < pp ) then
             ix = n - ix
           end if
           ignbin = ix
           return
         end if
   
         if ( ynorm + amaxp < alv ) then
           cycle
         end if
   
         x1 = real ( ix + 1, kind = 4 )
         f1 = fm + 1.0E+00
         z = real ( n + 1, kind = 4 ) - fm
         w = real ( n - ix + 1, kind = 4 )
         z2 = z * z
         x2 = x1 * x1
         f2 = f1 * f1
         w2 = w * w
   
         t = xm * log ( f1 / x1 ) + ( n - m + 0.5E+00 ) * log ( z / w ) &
           + real ( ix - m, kind = 4 ) * log ( w * p / ( x1 * q ) ) &
           + ( 13860.0E+00 - ( 462.0E+00 - ( 132.0E+00 - ( 99.0E+00 - 140.0E+00 &
           / f2 ) / f2 ) / f2 ) / f2 ) / f1 / 166320.0E+00 &
           + ( 13860.0E+00 - ( 462.0E+00 - ( 132.0E+00 - ( 99.0E+00 - 140.0E+00 &
           / z2 ) / z2 ) / z2 ) / z2 ) / z / 166320.0E+00 &
           + ( 13860.0E+00 - ( 462.0E+00 - ( 132.0E+00 - ( 99.0E+00 - 140.0E+00 &
           / x2 ) / x2 ) / x2 ) / x2 ) / x1 / 166320.0E+00 &
           + ( 13860.0E+00 - ( 462.0E+00 - ( 132.0E+00 - ( 99.0E+00 - 140.0E+00 &
           / w2 ) / w2 ) / w2 ) / w2 ) / w / 166320.0E+00
   
         if ( alv <= t ) then
           if ( 0.5E+00 < pp ) then
             ix = n - ix
           end if
           ignbin = ix
           return
         end if
   
       end if
   
     end do
   
     return
   end function

end module setting
