! To test random_number() performance with OpenMP

program main
use omp_lib
integer :: num_cpu = 1
real(8) :: t1,t2
real(4) :: r
integer :: i, j, len,status, nt = 10000
character*(8) :: str

!!$call get_command_argument(1, str, len, status)
!!$read(str,'(i2)') num_cpu
!!$write(*,*) 'Use num_cpu: ',num_cpu
!!$
!!$call omp_set_num_threads(num_cpu)

t1 = omp_get_wtime()

!$omp parallel do private(j,r)
do i = 1,10000
    do j = 1,nt
        call random_number(r)
        !print *, r
    enddo
enddo
!$omp end parallel do

t2 = omp_get_wtime()
write(*,*) 'Time: ',t2-t1
end

