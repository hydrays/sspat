module random
  implicit none

  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER, PARAMETER :: NPAR_CUMSUM=16
  INTEGER, PARAMETER :: NPAR_CUMPROD=8
  INTEGER, PARAMETER :: NPAR_POLY=8
  INTEGER, PARAMETER :: NPAR_POLYTERM=8
  REAL, PARAMETER :: PI = 3.1415927

  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
  INTEGER(K4B), SAVE :: lenran=0, seq=0
  INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
  INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
  INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
       nran,mran,ranv
  REAL, SAVE :: amm

  INTERFACE ran1
     module procedure ran1_s     
     module procedure ran1_v
  END INTERFACE

  INTERFACE ran2
     module procedure ran2_s     
     module procedure ran2_v
  END INTERFACE

  INTERFACE arth
     MODULE PROCEDURE arth_d, arth_i
  END INTERFACE

  INTERFACE gammln
     module procedure gammln_s
     module procedure gammln_v
  END INTERFACE

  INTERFACE expdev
     module procedure expdev_s
     module procedure expdev_v
  END INTERFACE

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm, &
          reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE

  INTERFACE ran_hash
     MODULE PROCEDURE ran_hash_s, ran_hash_v
  END INTERFACE

contains
  SUBROUTINE ran2_s(harvest)
    
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
    if (lenran < 1) call ran_init(1)
    rans=iran0-kran0
    if (rans < 0) rans=rans+2147483579_k4b
    iran0=jran0
    jran0=kran0
    kran0=rans
    nran0=ieor(nran0,ishft(nran0,13))
    nran0=ieor(nran0,ishft(nran0,-17))
    nran0=ieor(nran0,ishft(nran0,5))
    rans=iand(mran0,65535)
    mran0=ishft(3533*ishft(mran0,-16)+rans,16)+ &
         3533*rans+820265819_k4b
    rans=ieor(nran0,kran0)+mran0
    harvest=amm*merge(rans,not(rans), rans<0 )
  END SUBROUTINE ran2_s

  SUBROUTINE ran2_v(harvest)
    
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(OUT) :: harvest
    INTEGER(K4B) :: n
    n=size(harvest)
    if (lenran < n+1) call ran_init(n+1)
    ranv(1:n)=iran(1:n)-kran(1:n)
    where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
    iran(1:n)=jran(1:n)
    jran(1:n)=kran(1:n)
    kran(1:n)=ranv(1:n)
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
    ranv(1:n)=iand(mran(1:n),65535)
    mran(1:n)=ishft(3533*ishft(mran(1:n),-16)+ranv(1:n),16)+ &
         3533*ranv(1:n)+820265819_k4b
    ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)
    harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
  END SUBROUTINE ran2_v


  SUBROUTINE ran1_s(harvest)
    
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
    if (lenran < 1) call ran_init(1)
    rans=iran0-kran0
    if (rans < 0) rans=rans+2147483579_k4b
    iran0=jran0
    jran0=kran0
    kran0=rans
    nran0=ieor(nran0,ishft(nran0,13))
    nran0=ieor(nran0,ishft(nran0,-17))
    nran0=ieor(nran0,ishft(nran0,5))
    if (nran0 == 1) nran0=270369_k4b
    mran0=ieor(mran0,ishft(mran0,5))
    mran0=ieor(mran0,ishft(mran0,-13))
    mran0=ieor(mran0,ishft(mran0,6))
    rans=ieor(nran0,rans)+mran0
    harvest=amm*merge(rans,not(rans), rans<0 )
  END SUBROUTINE ran1_s

  SUBROUTINE ran1_v(harvest)
    
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(OUT) :: harvest
    INTEGER(K4B) :: n
    n=size(harvest)
    if (lenran < n+1) call ran_init(n+1)
    ranv(1:n)=iran(1:n)-kran(1:n)
    where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
    iran(1:n)=jran(1:n)
    jran(1:n)=kran(1:n)
    kran(1:n)=ranv(1:n)
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
    where (nran(1:n) == 1) nran(1:n)=270369_k4b
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
    ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
    harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
  END SUBROUTINE ran1_v

  SUBROUTINE gasdev(harvest)
    
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
    REAL :: rsq,v1,v2
    REAL, SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    if (gaus_stored) then
       harvest=g
       gaus_stored=.false.
    else
       do
          call ran1(v1)
          call ran1(v2)
          v1=2.0*v1-1.0
          v2=2.0*v2-1.0
          rsq=v1**2+v2**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq=sqrt(-2.0*log(rsq)/rsq)
       harvest=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev

  SUBROUTINE normdev(mean, var, harvest)
    
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
    REAL, INTENT(in) :: mean
    REAL, INTENT(in) :: var
    REAL :: g

    if (var < 0.0) then
       print *, 'error: var nagtive'
       read(*,*)
    else
       call gasdev(g)
       harvest=mean + var*g
    end if
  END SUBROUTINE normdev

  SUBROUTINE expdev_s(harvest)
    
    IMPLICIT NONE
    REAL, INTENT(OUT) :: harvest
    REAL :: dum
    call ran1(dum)
    harvest=-log(dum)
  END SUBROUTINE expdev_s

  SUBROUTINE expdev_v(harvest)
    
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(OUT) :: harvest
    REAL, DIMENSION(size(harvest)) :: dum
    call ran1(dum)
    harvest=-log(dum)
  END SUBROUTINE expdev_v

  FUNCTION poidev(xm)
    
    IMPLICIT NONE
    REAL, INTENT(IN) :: xm
    REAL :: poidev
    REAL :: em,harvest,t,y
    REAL, SAVE :: alxm,g,oldm=-1.0,sq
    if (xm < 12.0) then
       if (xm /= oldm) then
          oldm=xm
          g=exp(-xm)
       end if
       em=-1
       t=1.0
       do
          em=em+1.0
          call ran1(harvest)
          t=t*harvest
          if (t <= g) exit
       end do
    else
       if (xm /= oldm) then
          oldm=xm
          sq=sqrt(2.0*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.0)
       end if
       do
          do
             call ran1(harvest)
             y=tan(PI*harvest)
             em=sq*y+xm
             if (em >= 0.0)  exit
          end do
          em=int(em, 8)
          if (em+1.0 .le. 0.0) then
             print *, 'here'
             print *, xm
             print *, em
             read(*,*)
          end if
          t=0.9*(1.0+y**2)*exp(em*alxm-gammln_test(em+1.0)-g)
          call ran1(harvest)
          if (harvest <= t) exit
       end do
    end if
    poidev=em
  END FUNCTION poidev

  FUNCTION bnldev(pp,n)
    
    IMPLICIT NONE
    REAL, INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: n
    REAL :: bnldev
    INTEGER :: j
    INTEGER, SAVE :: nold=-1
    REAL :: am,em,g,h,p,sq,t,y,arr(24)
    REAL, SAVE :: pc,plog,pclog,en,oldg,pold=-1.0
    p=merge(pp,1.0-pp, pp <= 0.5 )
    am=n*p
    if (n < 25) then
       call ran1(arr(1:n))
       bnldev=count(arr(1:n)<p)
    else if (am < 1.0) then
       g=exp(-am)
       t=1.0
       do j=0,n
          call ran1(h)
          t=t*h
          if (t < g) exit
       end do
       bnldev=merge(j,n, j <= n)
    else
       if (n /= nold) then
          en=n
          oldg=gammln(en+1.0)
          nold=n
       end if
       if (p /= pold) then
          pc=1.0-p
          plog=log(p)
          pclog=log(pc)
          pold=p
       end if
       sq=sqrt(2.0*am*pc)
       do
          call ran1(h)
          y=tan(PI*h)
          em=sq*y+am
          if (em < 0.0 .or. em >= en+1.0) cycle
          em=int(em)
          t=1.2*sq*(1.0+y**2)*exp(oldg-gammln(em+1.0)-&
               gammln(en-em+1.0)+em*plog+(en-em)*pclog)
          call ran1(h)
          if (h <= t) exit
       end do
       bnldev=em
    end if
    if (p /= pp) bnldev=n-bnldev
  END FUNCTION bnldev

  FUNCTION arth_d(first,increment,n)
    REAL, INTENT(IN) :: first,increment
    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(n) :: arth_d
    INTEGER :: k,k2
    REAL :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER, INTENT(IN) :: first,increment,n
    INTEGER, DIMENSION(n) :: arth_i
    INTEGER :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  FUNCTION gammln_s(xx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: xx
    REAL :: gammln_s
    REAL :: tmp,x
    REAL :: stp = 2.5066282746310005
    REAL, DIMENSION(6) :: coef = (/76.18009172947146,&
         -86.50532032941677,24.01409824083091,&
         -1.231739572450155,0.1208650973866179e-2,&
         -0.5395239384953e-5/)
    if ( xx .le. 0.0 ) then
       print *, 'wrong argu of gamma'
       read(*,*)
    end if
    x=xx
    tmp=x+5.5
    tmp=(x+0.5)*log(tmp)-tmp
    gammln_s=tmp+log(stp*(1.000000000190015+&
         sum(coef(:)/arth(x+1.0,1.0,size(coef))))/x)
  END FUNCTION gammln_s

  FUNCTION gammln_test(xx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: xx
    REAL :: gammln_test
    REAL :: tmp,x
    REAL :: stp = 2.5066282746310005
    REAL, DIMENSION(6) :: coef = (/76.18009172947146,&
         -86.50532032941677,24.01409824083091,&
         -1.231739572450155,0.1208650973866179e-2,&
         -0.5395239384953e-5/)
    if ( xx .le. 0.0 ) then
       print *, 'test : wrong argu of gamma'
       read(*,*)
    end if
    x=xx
    tmp=x+5.5
    tmp=(x+0.5)*log(tmp)-tmp
    gammln_test=tmp+log(stp*(1.000000000190015+&
         sum(coef(:)/arth(x+1.0,1.0,size(coef))))/x)
  END FUNCTION gammln_test

  FUNCTION gammln_v(xx)
    IMPLICIT NONE
    INTEGER :: i
    REAL, DIMENSION(:), INTENT(IN) :: xx
    REAL, DIMENSION(size(xx)) :: gammln_v
    REAL, DIMENSION(size(xx)) :: ser,tmp,x,y
    REAL :: stp = 2.5066282746310005
    REAL, DIMENSION(6) :: coef = (/76.18009172947146,&
         -86.50532032941677,24.01409824083091,&
         -1.231739572450155,0.1208650973866179e-2,&
         -0.5395239384953e-5/)
    if (size(xx) == 0) RETURN
    if ( any(xx .le. 0.0) ) then
       print *, 'wrong argu of gamma'
       read(*,*)
    end if
    x=xx
    tmp=x+5.5
    tmp=(x+0.5)*log(tmp)-tmp
    ser=1.000000000190015
    y=x
    do i=1,size(coef)
       y=y+1.0
       ser=ser+coef(i)/y
    end do
    gammln_v=tmp+log(stp*ser/x)
  END FUNCTION gammln_v

  SUBROUTINE ran_init(length)
    IMPLICIT NONE
    INTEGER(K4B), INTENT(IN) :: length
    INTEGER(K4B) :: new,j,hgt
    if (length < lenran) RETURN
    hgt=hg
    if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
    if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
    if (hgt+1 /= hgng)    call nrerror('ran_init: arith assump 3 fails')
    if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
    if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
    if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
    if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
    if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
    if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
    if (lenran > 0) then
       ranseeds=>reallocate(ranseeds,length,5)
       ranv=>reallocate(ranv,length-1)
       new=lenran+1
    else
       allocate(ranseeds(length,5))
       allocate(ranv(length-1))
       new=1
       amm=nearest(1.0,-1.0)/hgng
       if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            call nrerror('ran_init: arth assump 10 fails')
    end if
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
    do j=1,4
       call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    end do
    where (ranseeds(new:,1:3) < 0) &
         ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
    where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
    if (new == 1) then
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
       rans=nran0
    end if
    if (length > 1) then
       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran
    end if
    lenran=length
  END SUBROUTINE ran_init
  !BL
  SUBROUTINE ran_deallocate
    if (lenran > 0) then
       deallocate(ranseeds,ranv)
       nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
       lenran = 0
    end if
  END SUBROUTINE ran_deallocate
  !BL
  SUBROUTINE ran_seed(sequence,size,put,get)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: sequence
    INTEGER, OPTIONAL, INTENT(OUT) :: size
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
    if (present(size)) then
       size=5*lenran
    else if (present(put)) then
       if (lenran == 0) RETURN
       ranseeds=reshape(put,shape(ranseeds))
       where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
       where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
    else if (present(get)) then
       if (lenran == 0) RETURN
       ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
       get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
       call ran_deallocate
       seq=sequence
    end if
  END SUBROUTINE ran_seed
  !BL
  SUBROUTINE ran_hash_s(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), INTENT(INOUT) :: il,ir
    INTEGER(K4B) :: is,j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_s
  !BL
  SUBROUTINE ran_hash_v(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
    INTEGER(K4B), DIMENSION(size(il)) :: is
    INTEGER(K4B) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_v

  FUNCTION reallocate_rv(p,n)
    REAL, DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv
  !BL
  FUNCTION reallocate_iv(p,n)
    INTEGER, DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv
  !BL
  FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
    INTEGER, INTENT(IN) :: n
    INTEGER :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_hv
  !BL
  FUNCTION reallocate_rm(p,n,m)
    REAL, DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER, INTENT(IN) :: n,m
    INTEGER :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_rm
  !BL
  FUNCTION reallocate_im(p,n,m)
    INTEGER, DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER, INTENT(IN) :: n,m
    INTEGER :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  SUBROUTINE assert(n, string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    logical, intent(in) :: n
    if (.not.(n)) then
       write (*,*) 'nrerror: asserteon failed ',string
       STOP 'program terminated... '
    end if
  END SUBROUTINE assert

end module random
