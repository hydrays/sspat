program ssa
  use nrtype
  use chem_data
  use random
  implicit none

  real(kind=8) x(NSpec)
  real(kind=8) te, u
  real(kind=8) tau, t, tp, td
  integer(I4B) is_nag
  integer(I4B) i, j, k, index

  call ran_seed(sequence=12341)
  !te = huge(1.0)
  te = 20000.0
  td = 5000.0
  do index = 1.0, NSample
     x = xinit
     t = 0.0
     tp = 0.0
     s1 = 1.0
     s2 = 1.0
     s3 = 0.005

     NP = 0.0
     NT = 0.0
     call getrate(x, a)
     do i = 1, NReac
        call expdev(u)
        NP(i) = u
     end do

     do while( .true. )

        if (t .ge. tp) then
           !call output_to_file(output_index)
           !call cell_stat(t)
           write(*, '((F12.4))', advance='no'), x(3)
           !output_index = output_index + 1
           tp = tp + 10.0
        end if

        if (t .ge. td) then
           !if ( mod(t-td, 1600.0) < 800 ) then
           !   s3 = 0.0
           !else
           !   s3 = 1.0
           !end if
           s1 = 0.0
           !s2 = 1.0
           !s3 = 1.0
        end if

        call Next_Reaction(k, tau)

        t = t + tau
        if ( t > te ) then
           !write(*, '(5(F12.4))'), t, x
           write(*,*)
           exit
        end if

        if ( k > 0 ) then
           x = x + nu(:, k)
           !if ( k .eq. 6 ) then
           !   u = 1.0
!              print *, 'here', k
           !else
              call expdev(u)
           !end if
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

!!$        if ( x(3).eq. 1.0 ) then
!!$           write(*, '(5(F12.4))'), t, x
!!$           exit
!!$        end if

     end do
  end do
end program ssa

subroutine checkx(x, is_nag)
  use chem_data
  use nrtype
  implicit none
  integer(I4B), intent(inout) :: x(NSpec)
  integer(I4B),  intent(out) :: is_nag
  is_nag = 0
  if(any(x < 0) ) then
     is_nag = 1
     print *, x
     print *, 'nag!'
  end if
end subroutine checkx

