program InvIsing
  implicit none
  integer, parameter :: p = 2 !/*system size*/
  real, parameter :: connprob = 1.0 !/*sparsity factor*/
  integer, parameter :: n = 10000	!/*number of samples*/
  integer, parameter :: JUMPS = 10000 !/*sampling period*/
  integer, parameter :: WARMUPS = 5000000 !/*MC warmups*/
  real, parameter :: B = 1.0 !/*inverse temperature*/
  real, parameter :: TOL = 0.000001 
  integer, parameter :: HES = 1

  integer INFO
  integer SEED
  integer i,j,k,l,r,ipiv(p-1),index(p),X(p)
  real C(p,p),m(p),Jtrue(p,p),Jest(p,p),Jnaive(p,p)
  real C(p,p),m(p),htrue(p),hest(p),hnaive(p)
  real delta
  real f(p-1),JAC(p-1,p-1)
  real gradnorm, gradmax
  real samples(n, p)
  integer samplev
  real thetas(p-1)
  integer done, count
  integer rnd
  real z, A1, A2, A3
  integer ii

  open(unit = 21, file="sample.dat", action = "read")
  do l = 1, n
     read(21, *) samples(l, :)
     !do k = 1, p
     !read(21, *) samplev
     !samples(l, k) = samplev
     !write(*, '(64(f4.0))'), samples(l, :)
     !end do
  end do
  close(21)
  !   printf("\n Done sampling \n");

  !   for(r=0;r<p;r++) {	
  do r = 1, p
     !print *, 'check 1', r
     done = 0
     count = 0
     do j = 1, p-1
        thetas(j) = 0.0
     end do
     do while (done .eq. 0)
        f = 0.0
        JAC = 0.0
        gradnorm = 0.0
        gradmax = 0.0
        !rnd= module(random(), HES)
        do l = 1, n
           z = 0.0
           do j = 1, p-1
              if (j .ge. r) then
                 ii = j + 1
              else
                 ii = j
              end if
              z = z + thetas(j)*samples(l, ii)
           end do
           z = exp(2.0*B*samples(l, r)*z)
           A1 = -1.0/n*2.0*B*samples(l, r)/(z+1.0)
           ! 	if ((l % HES) == rnd) {
           A2 = 1.0/n*(4.0*B*B*z)/((z+1.0)*(z+1.0))

           !print *, 'check 2', r

           do j = 1, p-1
              !print *, 'check 2.5', j
              if (j .ge. r) then
                 ii = j + 1
              else
                 ii = j
              end if
              f(j) = f(j) + samples(l, ii) * A1;	
              A3 = samples(l, ii)*A2;
              !print *, 'check 2.6', j			
              do k = 1, j
                 if (k .ge. r) then
                    ii = k + 1
                 else
                    ii = k
                 end if
                 !print *, 'check 2.61', j, k, l, ii			
                 JAC(j,k) = JAC(j,k) + samples(l, ii)*A3;
              end do
           end do
           ! 	}
           ! 	else {
           ! 	  for(j=0;j<p-1;j++) {
           ! 	    f[j]+=(((samples[l] >> j+(j>=r))&1)*2-1)*A1;	
           ! 	  }
           ! 	}
        end do

        !print *, 'check 3', r

        do j = 1, p-1
           do k = j+1, p-1
              ! 	for(k=j+1;k<p-1;k++) {
              JAC(j,k) = JAC(k,j);
           end do
           gradnorm = gradnorm + f(j)*f(j)
           if (gradmax < sqrt(f(j)*f(j))) gradmax=sqrt(f(j)*f(j))
        end do
        gradnorm = sqrt(gradnorm)

        write (*, "(2(f20.10))"), gradmax, gradnorm
        !       info = clapack_dgesv(CblasRowMajor, p-1, 1, &JAC[0][0], p-1, ipiv, &f[0], p-1);
        call DGESV(p-1, 1, JAC, p-1, ipiv, f, p-1, INFO)
        if (info .ne. 0) then
           write(*,*), "failure with error %d\n", info
           read(*,*)
        end if

        done = 1
        do j = 1, p-1
           thetas(j) = thetas(j) - f(j)
           if (f(j)>TOL .or. f(j)<-TOL) done = 0
        end do
        count = count + 1
     end do
     write(*,*) count, r
     do j = 1, p-1
        if (j .ge. r) then
           ii = j + 1
        else
           ii = j
        end if
        Jest(ii,r) = thetas(j)
     end do
     !     /*put thetas as columns in jest*/
  end do

  do j = 1, p
     Jest(j,j) = 0.0				
  end do



  !   /*calc delta*/
  !   for(j=0;j<p;j++) {
  !     for(i=0;i<j;i++) {
  !       delta+=((Jest[j][i]+Jest[i][j])/2.0-J[j][i])*((Jest[j][i]+Jest[i][j])/2.0-J[j][i]);
  !     }
  !   }
  !   delta=sqrt(delta/((p-1)/(2.0*connprob)));
  !   printf("%f\n",delta);



  !/*naive part*/
  m = 0.0
  C = 0.0
  do k=1, n
     do j=1, p
        m(j) = m(j) + samples(k, j)
        do i=1, p
           C(j,i) = C(j,i) + samples(k,j)*samples(k,i)
        end do
     end do
  end do
  do j=1, p
     m(j) = m(j)/n
  end do
  do j=1, p
     do i=1, p
        C(j,i) = C(j,i)/n
        C(j,i) = C(j,i) - m(j)*m(i)
     end do
  end do

  open(file="Bconnprob.dat",unit=31,action="write")
  write(31,"(2f10.2)"), B, connprob
  close(31)

  open(file="parametersest.dat",unit=101,action="write")
  do j = 1, p
     do i = 1, p
        write(101, "(f20.12)", advance="no"), Jest(j,i);		
     end do
     write(101,*)
  end do
  close(101);	

  ! open(file="parameters.dat",unit=102,action="write")
  ! do j = 1, p
  !    do i = 1, p
  !       write(102, "(f20.12)", advance="no"), Jtrue(j,i)
  !    end do
  !    write(102,*)
  ! end do
  ! close(102);	

  open(file="m.dat",unit=41,action="write")
  open(file="C.dat",unit=42,action="write")
  do j = 1, p
     write(41, "(f20.16)"), m(j)
     do i = 1, p
        write(42, "(f20.16)", advance="no"), C(j,i)
     end do
     write(42, *)
  end do
  close(41)
  close(42)

end program InvIsing
