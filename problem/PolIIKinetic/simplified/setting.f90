!!$ Pol II kinetic Simulation

!!$ >>> NSpec = 3

!!$ (1,  2,  3)
!!$  C   R   Vita

!!$ >>> NReac = 6

!!$  1: I -> C
!!$  2: C -> I
!!$  3: C -> R
!!$  4: R -> C
!!$  5: R -> E + C
!!$  6: R -> E

module setting
  implicit none
  integer, parameter :: NSample = 2000
  integer, parameter :: NSpec=3
  integer, parameter :: NReac=6
  real, parameter :: ep = 0.1
  integer, parameter :: E0 = 2
  integer, parameter :: E1 = 10
  real, parameter :: v_mature = 100.0
!  real, parameter :: t_on = 1460.0
!  real, parameter :: t_off = 1460.0

!  real, parameter :: t_on = 0.5
!  real, parameter :: t_off = 0.5

  real, parameter :: t_on = 840000.0
  real, parameter :: t_off = 84.0

  real env
  real s(3)
  real a(NReac)
  real NT(NReac)
  real NP(NReac)
  real pop_ratio(3)
  real env_t

  integer, parameter, dimension(NSpec,NReac) :: nu = reshape( &
       source=(/ &
       (/01, 00, -E0/), & !1 
       (/-1, 00, 00/), & !2
       (/-1, 01, -E0/), & !3
       (/01, -1, 00/), & !4
       (/01, -1, -E0/), & !5
       (/00, -1, -E0/) & !6
       /), shape = (/NSpec, NReac/) &
       )

  real, parameter, dimension(NReac) :: c =  &
       (/0.1, 0.02, 1.0, 0.2, 10.0*(1.0-ep), 10.0*ep/)
       !(/0.1, 0.05, 1.0, 0.1, 10.0*(1.0-ep), 10.0*ep/)

  real, parameter, dimension(NSpec) :: xinit =  &
       (/0.0, 0.0, 50.0/)

  type cell
     integer id
     real x(NSpec)
  end type cell

  type(cell) CellPool(2000)

contains
  subroutine getrate(x, a)
    implicit none
    real, intent(in) :: x(NSpec)
    real, intent(out) :: a(NReac)
    real is_avai, act
    if ( (x(1) + x(2)) .eq. 1.0) then
       is_avai = 0.0
    else if ( (x(1) + x(2)) .eq. 0.0) then
       is_avai = 1.0
    else
       is_avai = -1.0
       print *, 'wrong occupy'
       print *, x
       read(*,*)
    end if
    
    if ( x(3) .le. E0 ) then
       act = 0.0
    else
       act = 1.0
    end if

    a(1) = act*s(1)*c(1)*is_avai
    a(2) = c(2)*x(1)
    a(3) = act*s(2)*c(3)*x(1)
    a(4) = c(4)*x(2)
    a(5) = act*s(3)*c(5)*x(2)
    a(6) = act*s(3)*c(6)*x(2)
  end subroutine getrate

  subroutine Next_Reaction(k, tau)
    implicit none
    integer, intent(out) :: k
    real, intent(out) :: tau
    real tau_temp
    integer i

    tau = huge(0.0)
    k = 0
    do i = 1, NReac
       if ( a(i) > 0.0 ) then
          tau_temp = ( NP(i) - NT(i) ) / a(i)
          if ( tau_temp < tau) then
             tau = tau_temp
             k = i
          end if
       else
       end if
    end do
    if ( tau > 0.1 ) then
       tau = 0.1
       k = 0
    end if
    !if ( k .eq. 0 ) then
    !   tau = 0.1
    !end if
    if ( tau < 0 .or. k < 0 ) then
       write(*,*), 'error', 'tau', tau
       print *, k, NP(k) - NT(k), a(k)
       read(*,*)
    end if
  end subroutine Next_Reaction

  subroutine evolve_cell(index, te)
    use random
    implicit none
    integer, intent(in) :: index
    real, intent(in) :: te
    real t, u, tau
    integer i, j, k
    real x(NSpec)
    real cs

    x = CellPool(index)%x
    if ( env > 0.5 ) then 
       cs = 1.0
    else
       cs = 0.0
    end if
    if ( CellPool(index)%id .eq. 3 ) then
       s = (/1.0, 1.0, cs/)
    else if ( CellPool(index)%id .eq. 2 ) then
       s = (/1.0, cs, 1.0/)
    else if ( CellPool(index)%id .eq. 1 ) then
       s = (/cs, 1.0, 1.0/)
    else
       print *, "error in choosing s... stop"
       read(*,*)
       stop
    end if

    t = 0.0
    NP = 0.0
    NT = 0.0
    call getrate(x, a)
    do i = 1, NReac
       call expdev(u)
       NP(i) = u
    end do
    do while( t < te )
       call Next_Reaction(k, tau)
       t = t + tau
       if ( t > te ) then
          CellPool(index)%x = x
          exit
       end if
       if ( k > 0 ) then
          x = x + nu(:, k)
          ! For model that consider Vita (x(5)) only. 
          if ( k.eq.5 .or. k.eq.6 ) then
             x(3) = x(3) + env * real(E1)
          end if
          call expdev(u)
          if ( u.eq.0 ) then
             print *, 'exp value = 0'
             read(*,*)
          end if
          NP(k) = NP(k) + u
       end if
       do i = 1, NReac
          NT(i) = NT(i) + a(i)*tau
       end do
       call getrate(x, a)
    end do
  end subroutine evolve_cell

  subroutine output_to_file(index)
    implicit none

    integer, intent(in) :: index
    character(30) filename
    integer i

    WRITE(filename,'(A7,I5.5,A4)') './out/m', index, '.dat'
    open (unit = 11, file=filename, action="write")

    do i = 1, NSample
       write(11, '(I5, 10(f10.2))'), CellPool(i)%id, CellPool(i)%x
    end do
    close(11)
  end subroutine output_to_file

  subroutine output_to_file2(unit, time)
    implicit none
    integer, intent(in) :: unit
    real, intent(in) :: time
    integer i
    integer np

    write(unit, '(f10.2, f10.2, 3(I5, f10.2), 3(f10.2))'), time, env, &
         CellPool(1)%id, CellPool(1)%x(3), &
         CellPool(2)%id, CellPool(2)%x(3), &
         CellPool(3)%id, CellPool(3)%x(3), &
         pop_ratio
    
  end subroutine output_to_file2

  subroutine checkx(x, is_nag)
    implicit none
    integer, intent(inout) :: x(NSpec)
    integer,  intent(out) :: is_nag
    is_nag = 0
    if(any(x < 0) ) then
       is_nag = 1
       print *, x
       print *, 'nag!'
    end if
  end subroutine checkx

  subroutine update_env(t)
    use random
    implicit none
    real, intent(in) :: t
    real u
    call gasdev(u)
    ! update scheme 1
    !env = 10.0
    
    ! update scheme 2
    if (t .eq. 0.0) then
       env_t = 0.0
    end if

    if ( t .ge. env_t ) then
       if ( env .eq. 0.0 ) then
          env = 1.0
          env_t = env_t + t_on
       else if ( env .eq. 1.0 ) then
          env = 0.0
          env_t = env_t + t_off
       else
          print *, "something wrong"
          read(*,*)
       end if
    end if
    ! update scheme 3
    !env = max(1.0, 5.0*(1.0+cos(0.001*t)) + 0.5*u)
  end subroutine update_env

end module setting
