
      PROGRAM  KMC_DIST_TIMING

       use main
       IMPLICIT NONE

       INTEGER(KIND=int16) :: i,j,i1,jj,j1,Mx,My,ncell,err,jmax,n,n21,n2,ib,i2,i3,ix,iy,iix,iiy,nruns,nbins
       INTEGER(KIND=int16) :: ix1,iy1,nlvls,nlvlg,natmsw,nSiatms,nSiatms2,nGeatms,nsfatms,nomove,nomove1,frz1y,sfid,seedrng
       INTEGER(KIND=int16) :: ii,pp,nn,mm,tt,ss,numf,nsf1,numnbf,ndiv,naccept,ndepst,sdsteps,count,numout,startflag,numstart
       REAL(KIND=r8) :: xix,yix,dxs2,dys,sigs,sigg,epss,epsg,boxl2,boxlen
       REAL(KIND=r8) :: yt,tol,mina,maxa,rr,dxd,dxsrad,dxd_br,MinEne,PrevMinEne,Uene,Uene1,Uene2,Uene3
       REAL(KIND=r8) :: ExRate,ratio,yy,Enediff,MineneT,nbf,deprate,depflux,slope,TempK
       REAL(KIND=r8) :: tempsig,tempeps,misfit,radnc,angle,dtheta,maxheight,xitemp,yitemp,xitempn,yitempn
       REAL(KIND=r8) :: beta,factor,SimulTime1,SimulTime2,Udist,ratetot,rt,ymax,rnum,st0,st1,st2,st3,st4
       REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xx,xxtemp
       REAL(KIND=r8),PARAMETER :: dEAdjust=1.411050034856138_r8
       REAL(KIND=r8),PARAMETER :: dxs=3.0194_r8 !6.0956_r8
       REAL(KIND=r8),PARAMETER :: kB=8.617342311969807e-05_r8
       REAL(KIND=r8),PARAMETER :: kfac=1.0e12_r8
       REAL(KIND=r8),PARAMETER :: intercept=0.0_r8 
       REAL(KIND=r8),PARAMETER :: twopi=6.28318530717959_r8
       INTEGER(KIND=int16),PARAMETER :: radatm=4_int16
       INTEGER(KIND=int16),DIMENSION(1,2) :: Abin
       INTEGER(KIND=int16),ALLOCATABLE,DIMENSION(:,:) :: Gbin,cdnumtemp
       REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: Varrt,sumrates,Strain,Rates,RatesTemp
       REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: Strainxy
       INTEGER(KIND=int16),ALLOCATABLE,DIMENSION(:) :: Gbin2,MinLocEne
       INTEGER(KIND=int16),DIMENSION(8) :: time1,time2,time3,time4,time5,time6,time7,time8,time9,time10,time11,time12 
       CHARACTER (LEN=20) :: sampin
       CHARACTER (LEN=8)  :: datec
       CHARACTER (LEN=10) :: timec
       CHARACTER (LEN=5)  :: zonec
       LOGICAL :: intst,biehl 
       REAL(KIND=r8)      ::totaltime,subtime1
       REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) ::jac
       REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: eigvec
       REAL(KIND=r8)::eigval
       REAL(KIND=r8), allocatable, dimension(:,:) ::lower, upper
       Integer:: INFO
       double precision, allocatable, dimension(:,:):: VR
       double precision, allocatable, dimension(:):: WR,WI,WORK,VL

       READ(5,*)startflag,numout,Mx,My,seedrng,factor,misfit
       READ(5,*)natmsw,nlvls,nlvlg,TempK
       READ(5,*)sampin
       READ(5,*)natms,nGeatms,numstart,naccept
       READ(5,*)boxlen, depflux, tol, nruns

! Initialize variables
       beta=1.0_r8/(kB*TempK)
       frzatms1=natmsw+1_int16
       epss=0.4_r8; epsg=0.3387_r8
       slope=bounds(misfit)
       deprate=depflux*REAL(natmsw,KIND=r8)
       sigs=2.7153_r8; sigg=sigs*(1.0_r8+misfit); FmaxTol=1.0e-02_r8
       IF(misfit.EQ.0.0_r8)epsg=epss
          
       dxs2=0.5_r8*dxs; dys=0.5_r8*dxs*SQRT(3.0_r8); dys2=0.5_r8*dys
       ndiv=6_int16; frz1y=natms+frzatms1; ss=10000_int16; numnbf=0_int16
       latdist2=2.0_r8*sigg
       !cutoff radius
       rcut=4.0_r8*dxs; rcut2=rcut*rcut; rrcut2=1.0_r8/rcut2; rrcut6=rrcut2*rrcut2*rrcut2
       dxss=dxs !Globalize dxs
       dxd=1.5_r8*dys !sqrt(3.0_r8)*dxs-0.1_r8
       dxsrad=REAL(radatm,KIND=r8)*dxs

       epss=factor*epss; epsg=factor*epsg; FmaxTol=factor*FmaxTol

       nSiatms=natmsw*nlvls;nSiatms2=2_int16*nSiatms 
       natms2=2_int16*natms; jmax=10_int16; ncell=Mx*My; mapsiz=8_int16*ncell 

       ALLOCATE(xx(natms2),xxtemp(natms2),fp(natms2),f2p(natms2,natms2),NewR(ndiv,2),&
                sig(natms),eps(natms),ngbharr(natms),ngbharr_br(natms),ngbh1arr_br(natms),map(mapsiz),&
                HEAD(ncell),HEADT(ncell),LIST(natms),LISTT(natms),COORDNUM(natms,2),&
                cdnumtemp(natms,2),EneArr(ndiv),MinLocEne(ndiv),Strain(natms),Strainxy(natms,2),&
                RatesTemp(natms),Rates(natms),jac(natms2,natms2),eigvec(natms2),lower(natms2,natms2),upper(natms2,natms2),&
                VL(natms2),WR(natms2),WI(natms2),VR(natms2,natms2),WORK(4*natms2),  STAT=err)
         IF(err .NE. 0) THEN
           WRITE(6,*)'Error allocating data arrays :1'
         END IF


!coordinates for generating new positions
       NewR(1,1)=dxs;           NewR(1,2)=0.0_r8
       NewR(2,1)=0.5_r8*dxs;    NewR(2,2)=dys
       NewR(3,1)=0.5_r8*dxs;    NewR(3,2)=-1.0_r8*dys
       NewR(4,1)=-0.5_r8*dxs;   NewR(4,2)=dys
       NewR(5,1)=-1.0_r8*dxs;   NewR(5,2)=0.0_r8
       NewR(6,1)=-0.5_r8*dxs;   NewR(6,2)=-1.0_r8*dys

       CALL SETUP_MAP(Mx,My)

   IF(startflag==0) THEN
! Set initial boxlength. This is just for creating initial configuration. It will be updated
! later
     boxl=REAL(natmsw-1,KIND=r8)*dxs
! Generate initial configuration and store in xx. Also store lattice distance and well-depth
   ! Silicon substrate
     i1=0
     DO j=1,nlvls
        j1=j-1 
        yt=REAL(j1,KIND=r8)*dys
        boxl2=-0.5_r8*boxl+dxs2*REAL(MOD(j1,2),KIND=r8)
        DO i=1,natmsw
           i1=i1+1; i2=i1+natms
           xx(i1)=boxl2+REAL(i-1,KIND=r8)*dxs; xx(i2)=yt
           sig(i1)=sigs; eps(i1)=epss
        END DO
      END DO

   ! Germanium layer
     DO j=1,nlvlg
        j1=nlvls+j-1
        yt=REAL(j1,KIND=r8)*dys
        boxl2=-0.5_r8*boxl+dxs2*REAL(MOD(j1,2),KIND=r8)
        DO i=1,natmsw
           i1=i1+1; i2=i1+natms
           xx(i1)=boxl2+REAL(i-1,KIND=r8)*dxs; xx(i2)=yt
           sig(i1)=sigg; eps(i1)=epsg
        END DO
      END DO

!  Write initial configuration into file

       OPEN(unit=10,file='init_conf.txt',status='unknown',action='write')
       DO i=1,natms2
          WRITE(10,16)xx(i)
       END DO

       CLOSE(10)

        boxl=boxlen
        boxl2=0.5_r8*boxl
        rboxl=1.0_r8/boxl
	
	CALL HESSIAN(Uene,xx,jac)
	open(unit=26, file='HESSIAN.txt', ACTION="write", STATUS="replace")
	write(26,*) int( natms2 )
	do i=1,natms2
  	    write(26,*) real( jac(i,:) ) 
	end do
	close(26)


!	CALL DominantEigen(jac,eigvec,eigval)
!	print*,'initial dominant eigenvalue:', eigval
!	Call LUfactor(jac,lower,upper)

!	open(unit=28, file='Lower.txt', ACTION="write", STATUS="replace")
!        do i=1,natms2
!            write(28,*) real( lower(i,:) )
!        end do
!	close(28)
	
!	open(unit=29, file='Upper.txt', ACTION="write", STATUS="replace")
!	do i=1,natms2
!            write(29,*) real( upper(i,:) )
!        end do
!        close(29)

!	Call LeastDominant (lower,upper,eigvec,eigval)
!	print*,'initial weak eigenvalue:', eigval

!	Call DGEEV('N','V',natms2,jac,natms2,WR,WI,VL,natms2,VR,natms2,WORK,4*natms2,INFO)
!	print*, WR
        CALL CONJ_GRAD(xx,tol,jmax,Enediff)
        CALL HESSIAN(Uene,xx,jac)

	open(unit=27, file='RelaxedHESSIAN.txt', ACTION="write", STATUS="replace")
	write(27,*) int( natms2 )
        do i=1,natms2
            write(27,*) real( jac(i,:) )
        end do
        close(27)


!	CALL DominantEigen(jac,eigvec,eigval)
!	print*,'relaxed dominant eigenvalue', eigval
!	Call LUfactor(jac,lower,upper)
!	Call LeastDominant (lower,upper,eigvec,eigval)	
!	print*,'relaxed weak eigenvalue:', eigval

	CALL FULL_POT_ENE(xx,MinEne)
        CALL LINKS(Mx,My,xx)

        !boxl=maxval(xx(1:natms))-minval(xx(1:natms)) + dxs2 
        !boxl2=0.5_r8*boxl;
        !rboxl=1.0_r8/boxl;

        WRITE(6,*)maxval(xx(1:natms))-minval(xx(1:natms)) + dxs2
        dxd_br=dxd+4.0_r8*dxs
        DO i=1,natms

           Abin(1,1)=i; Abin(1,2)=0_int16
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
           CALL POT_DIST(i,XX,Udist,Uene)
           Rates(i)=kfac*exp(beta*(slope*Udist+intercept+Uene))

        END DO



        OPEN(unit=8,file='relaxed_iniconf.txt',status='unknown',action='write')
        DO i=1,natms2
           WRITE(8,16)xx(i)
        END DO

        CLOSE(8)

!Find coordination number for each atom
        DO i=1,natms
           COORDNUM(i,1)=i
        END DO

        COORDNUM(:,2)=0_int16
        CALL LINKS(Mx,My,xx)
        !dxd=dxs; 
        dxd_br=dxd+3.0_r8*dxs
        CALL COORD_NUM(natms,COORDNUM,Mx,My,dxd,dxd_br,xx)
        
        ndepst=0_int16
   ELSE

       boxl=boxlen 
       boxl2=0.5_r8*boxl;
       rboxl=1.0_r8/boxl;
   
       ndepst=nGeatms
       OPEN(unit=9,file=sampin,status='old',action='read')
       DO i=1,natms2
          READ(9,*)xx(i)
       END DO

       CLOSE(9)

       DO i=1,nSiatms
          sig(i)=sigs; eps(i)=epss
       END DO

       DO i=nSiatms+1,natms
          sig(i)=sigg; eps(i)=epsg
       END DO

!Global relaxation
        CALL CONJ_GRAD(xx,tol,jmax,Enediff)
        CALL FULL_POT_ENE(xx,MinEne)
        CALL LINKS(Mx,My,xx)

        dxd_br=dxd+4.0_r8*dxs
        DO i=1,natms

           Abin(1,1)=i; Abin(1,2)=0_int16
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
           CALL POT_DIST(i,XX,Udist,Uene)
           Rates(i)=kfac*exp(beta*(slope*Udist+intercept+Uene))

        END DO

        OPEN(unit=8,file='relaxed_iniconf.txt',status='unknown',action='write')
        DO i=1,natms2
           WRITE(8,16)xx(i)
        END DO

        CLOSE(8)

!Find coordination number for each atom
        DO i=1,natms
           COORDNUM(i,1)=i
        END DO

        COORDNUM(:,2)=0_int16
        CALL LINKS(Mx,My,xx)
        !dxd=dxs; 
        dxd_br=dxd+3.0_r8*dxs
        CALL COORD_NUM(natms,COORDNUM,Mx,My,dxd,dxd_br,xx)

   END IF


      n=natms2; n2=natms; n21=n2+1

       nomove=0_int16; nomove1=0_int16
  SimulTime1=0.0_r8; SimulTime2=0.0_r8    
  CALL RANDOM_SEED 
  DO i=1,seedrng
     CALL RANDOM_NUMBER(rnum)
  END DO

  ii=numout; pp=numout+10000
  numf=numout+1000_int16
  DO jj=numstart,nruns
     numsteps=jj 
! Perform a global relaxation of the system if i=1 or i=1000000

  CALL DATE_AND_TIME(datec,timec,zonec,time1)

      IF(MOD(jj,100000_int16)==0_int16)THEN
         CALL CONJ_GRAD(xx,tol,jmax,Enediff)
         CALL FULL_POT_ENE(xx,MinEne)

         OPEN(unit=10,file='enecheck.txt',status='UNKNOWN',position='APPEND',action='write')
         WRITE(10,'(I9,2X,F20.12,2X,F20.12)')jj,MinEne,MinEne/REAL(natms,KIND=r8)
         close(10)

         COORDNUM(:,2)=0_int16
         CALL LINKS(Mx,My,xx)
         !dxd=dxs; 
         dxd_br=dxd+4.0_r8*dxs
         CALL COORD_NUM(natms,COORDNUM,Mx,My,dxd,dxd_br,xx)
         DO i=1,natms
            Abin(1,1)=i; Abin(1,2)=0_int16
            CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
            CALL POT_DIST(i,XX,Udist,Uene)
            Rates(i)=kfac*exp(beta*(slope*Udist+intercept+Uene))
         END DO

! Write relaxed configuration to a file and compute strain
         ii=ii+1_int16
         OPEN(unit=ii,status='UNKNOWN',action='write')
         WRITE(ii,16)xx
         close(ii)
      END IF

      xxtemp=xx
      cdnumtemp=COORDNUM
      PrevMinEne=MinEne
 
! Find the surface atoms

      mina=minval(xx(1:natms))
      maxa=maxval(xx(1:natms))
    
     
      nsfatms=0_int16
      DO i=frzatms1,natms
         IF(COORDNUM(i,2) < 5)nsfatms=nsfatms+1
      END DO

      ALLOCATE(Varrt(nsfatms),Gbin2(nsfatms),sumrates(nsfatms),STAT=err)
      IF(err .NE. 0) THEN
        WRITE(6,*)'Error allocating Varrt, Gbin2 and sumrates arrays'
      END IF

      nn=0_int16
      DO i=frzatms1,natms
         IF(COORDNUM(i,2) < 5)THEN
            nn=nn+1_int16
            Gbin2(nn)=COORDNUM(i,1)
         END IF
      END DO

      Varrt=0.0_r8         
      !dxd=dxs; dxd_br=dxd+4.0_r8*dxs
      CALL LINKS(Mx,My,XX)

      IF(MOD(jj,100000_int16)==0_int16)THEN
         pp=pp+1_int16
         OPEN(unit=pp,status='UNKNOWN',action='write')
      END IF

      DO i=1,nsfatms
         ix=Gbin2(i)
         Varrt(i)=Rates(ix)
         IF(MOD(jj,100000_int16)==0_int16)WRITE(pp,18)i,xx(ix),xx(ix+natms)
      END DO

      close(pp)

     DO i=1,nsfatms
       sumrates(i)=sum(Varrt(1:i))
     END DO

!Sum up approximate hopping rates
     ratetot=sumrates(nsfatms) + deprate !depflux*REAL(nsfatms,KIND=r8)

     CALL RANDOM_NUMBER(rnum)
     rt=ratetot*rnum

! Maximum y coordinate
     ymax=MAXVAL(XX(n21:n))
     maxheight=ymax+1.2_r8*dys
      
! Search sumrates to find out which surface atom to move or whether to deposit a Ge
! atom
     
     CALL BINARY_SEARCH(rt,sumrates,sfid)
     IF(sfid .NE. 0) THEN !We will attempt to move the corresponding atom

       ix=Gbin2(sfid); iy=ix+natms

!Update the rates of the atoms around the atom to be moved. That is eliminate
!the effect of the atom to be moved on the neighboring atoms

        tempsig=sig(ix); tempeps=eps(ix)
        sig(ix)=0.0_r8; eps(ix)=0.0_r8
        dxd_br=dxd+4.0_r8*dxs
        Abin(1,1)=ix; Abin(1,2)=0_int16
        CALL LINKS(Mx,My,XX)
        CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)

        intbrnatms=intnatms+brnatms
        ALLOCATE(Gbin(intbrnatms,2),STAT=err)
        IF(err .NE. 0) THEN
          WRITE(6,*)'Error allocating Gbin array'
        END IF

        nn=0_int16
        DO i=1,intnatms
           nn=nn+1_int16
           Gbin(nn,1)=ngbharr(i); Gbin(nn,2)=0_int16
        END DO

        DO i=1,brnatms
           nn=nn+1_int16
           Gbin(nn,1)=ngbharr_br(i); Gbin(nn,2)=0_int16
        END DO

        DO i=1,intbrnatms
           nn=Gbin(i,1)
!update rates
           Abin(1,1)=nn; Abin(1,2)=0_int16
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
           CALL POT_DIST(nn,XX,Udist,Uene)
           Rates(nn)=kfac*exp(beta*(slope*Udist+intercept+Uene))
        END DO
        DEALLOCATE(Gbin)

        sig(ix)=tempsig; eps(ix)=tempeps

!Find atoms within 4 lattice spaces of atom ix

        dxd_br=dxd+4.0_r8*dxs
        Abin(1,1)=ix; Abin(1,2)=0_int16
        CALL LINKS(Mx,My,XX)
        CALL COORD_NUM(1_int16,Abin,Mx,My,4.0*dxs,dxd_br,xx)
        ALLOCATE(Gbin(intnatms,2),STAT=err)
        IF(err .NE. 0) THEN
          WRITE(6,*)'Error allocating Gbin array'
        END IF

        nn=0_int16
        DO i=1,intnatms
           IF(COORDNUM(ngbharr(i),2) .LE. 4_int16) THEN
              nn=nn+1_int16
              Gbin(nn,1)=ngbharr(i)
           END IF 
        END DO
        IF(nn > 0_int16) THEN
           naccept=naccept+1_int16
           mm=1_int16+NINT(rnum*REAL(nn-1_int16),KIND=int16)
           iix=Gbin(mm,1); iiy=iix+natms
           xitemp=XX(ix); yitemp=XX(iy)
           CALL NEW_COORD(ix,iix,ndiv,Mx,My,dxs,xx)
         END IF
         DEALLOCATE(Gbin)

           IF(MINVAL(EneArr) < 0.0_r8) THEN
              MinLocEne=MinLoc(EneArr); mm=MinLocEne(1)
              XX(ix)=XX(iix)+NewR(mm,1); XX(iy)=XX(iiy)+NewR(mm,2)
              XX(ix)=XX(ix)-boxl*ANINT(XX(ix)*rboxl)
           ELSE
              OPEN(unit=2,file='no_hop.txt',status='UNKNOWN',position='REWIND',action='write')
              nomove=nomove+1_int16
              WRITE(2,*)nomove,naccept,jj
              WRITE(2,*)ix,xitemp,yitemp
              WRITE(2,16)xx
              close(2)

           END IF
           
           dxd_br=2.0_r8*dxs
           xitempn=XX(ix); yitempn=XX(iy)
           Abin(1,1)=ix; Abin(1,2)=0_int16
           !dxd=REAL(radatm,KIND=r8)*dxs; dxd_br=dxd+4.0_r8*dxs
           CALL LINKS(Mx,My,XX)
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxs,dxd_br,xx)
           DO WHILE (brnatms==0)
              xx(iy)=xx(iy)-dys
              CALL COORD_NUM(1_int16,Abin,Mx,My,dxs,dxd_br,xx)
           END DO

           CALL STEEP_DESC(ix,xx,tol,sdsteps)

           dxd_br=dxsrad+4.0_r8*dxs
           Abin(1,2)=0_int16
           CALL LINKS(Mx,My,XX)
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxsrad,dxd_br,xx)
           CALL CONJ_GRAD_SUBSET(2*intnatms,jmax,xx,tol)
           CALL POT_1_NGBH(XX,Uene2)
           MinEne=MinEne-Uene1+Uene2
 
           CALL BOUNDARY_FORCE(XX,nbf,intst)

           IF(MinEne > 0.0_r8)THEN
              CALL FULL_POT_ENE(xx,MinEne)
           END IF

           IF( MinEne > 0.0_r8 ) THEN
              xx=xxtemp
              CALL CONJ_GRAD(xx,tol,jmax,Enediff)
              CALL FULL_POT_ENE(xx,MinEne)

              COORDNUM(:,2)=0_int16
              CALL LINKS(Mx,My,xx)
              !dxd=dxs; dxd_br=dxd+4.0_r8*dxs
              CALL COORD_NUM(natms,COORDNUM,Mx,My,dxd,dxd_br,xx)
              DO i=1,natms

                 Abin(1,1)=i; Abin(1,2)=0_int16
                 CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
                 CALL POT_DIST(i,XX,Udist,Uene)
                 Rates(i)=kfac*exp(beta*(slope*Udist+intercept+Uene))

              END DO

           END IF
      
     ELSE ! Deposition

        ndepst=ndepst+1 
        nGeatms=nGeatms+1; natms=natms+1; natms2=2*natms
        frz1y=natms+frzatms1
        RatesTemp=Rates
 
        DEALLOCATE(xx,sig,eps,fp,f2p,ngbharr,ngbharr_br,ngbh1arr_br,LIST,LISTT,COORDNUM,Strain,Strainxy,Rates)
       
        ALLOCATE(xx(natms2),fp(natms2),f2p(natms2,natms2),sig(natms),&
                 eps(natms),ngbharr(natms),ngbharr_br(natms),ngbh1arr_br(natms),LIST(natms),&
                 LISTT(natms),COORDNUM(natms,2),Strain(natms),Strainxy(natms,2),Rates(natms),STAT=err)
        IF(err .NE. 0) THEN
           WRITE(6,*)'Error allocating data arrays :2'
        END IF

        sig(1:nSiatms)=sigs; sig(nSiatms+1:natms)=sigg
        eps(1:nSiatms)=epss; eps(nSiatms+1:natms)=epsg
 
        xx(1:n2)=xxtemp(1:n2)
        xx(natms+1:natms+n2)=xxtemp(n21:n)
        COORDNUM(1:n2,1)=cdnumtemp(1:n2,1)
        COORDNUM(1:n2,2)=cdnumtemp(1:n2,2)
        COORDNUM(natms,1)=natms
        COORDNUM(natms,2)=0_int16
        Rates(1:n2)=RatesTemp; Rates(natms)=0.0_r8
        ix=natms; iy=ix+natms

        !dxd=2.0_r8*REAL(radatm,KIND=r8)*dxs; dxd_br=dxd+4.0_r8*dxs
        dxd_br=2.0_r8*dxs
        Abin(1,1)=ix 

           CALL RANDOM_NUMBER(rnum)
           xix=mina+(maxa-mina)*rnum 
           xix=xix-boxl*ANINT(xix*rboxl)
           yix=ymax+dxs
           xx(ix)=xix; xx(iy)=yix
           Abin(1,2)=0_int16
           CALL COORD_NUM_SIMPLE(1_int16,Abin,dxs,dxd_br,xx)
           DO WHILE (brnatms==0)
              xx(iy)=xx(iy)-dys
              CALL COORD_NUM_SIMPLE(1_int16,Abin,dxs,dxd_br,xx)
           END DO
           CALL STEEP_DESC(ix,xx,tol,sdsteps)

        !dxd=REAL(radatm,KIND=r8)*dxs; dxd_br=dxd+4.0_r8*dxs
           dxd_br=dxsrad+4.0_r8
           Abin(1,2)=0_int16
           CALL COORD_NUM_SIMPLE(1_int16,Abin,dxsrad,dxd_br,xx)
           CALL DEPOSIT_LOCAL_RELAX(2*intnatms,jmax,xx,tol)
           CALL POT_DIST(ix,XX,Udist,Uene3)
           WRITE(9,*)ndepst,Uene3

          n=natms2; n2=natms; n21=n2+1

          CALL BOUNDARY_FORCE(XX,nbf,intst)


        DEALLOCATE(xxtemp,cdnumtemp,RatesTemp)

        ALLOCATE(xxtemp(natms2),cdnumtemp(natms,2),RatesTemp(natms),STAT=err)
        IF(err .NE. 0) THEN
           WRITE(6,*)'Error allocating xxtemp'
        END IF

     END IF

     DEALLOCATE(Varrt,sumrates)

!Update COORDNUM
        intbrnatms=intnatms+brnatms
        ALLOCATE(Gbin(intbrnatms,2),STAT=err)
        IF(err .NE. 0) THEN
          WRITE(6,*)'Error allocating Gbin array'
        END IF

        nn=0_int16
        DO i=1,intnatms
           nn=nn+1_int16
           Gbin(nn,1)=ngbharr(i); Gbin(nn,2)=0_int16
        END DO

        DO i=1,brnatms
           nn=nn+1_int16
           Gbin(nn,1)=ngbharr_br(i); Gbin(nn,2)=0_int16
        END DO
        CALL LINKS(Mx,My,xx)
        CALL COORD_NUM(intbrnatms,Gbin,Mx,My,dxd,dxd_br,xx)
        dxd_br=dxd+4.0_r8*dxs
        DO i=1,intbrnatms
           nn=Gbin(i,1)
           COORDNUM(nn,2)=Gbin(i,2)
!update rates
           Abin(1,1)=nn; Abin(1,2)=0_int16
           CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
           CALL POT_DIST(nn,XX,Udist,Uene)
           Rates(nn)=kfac*exp(beta*(slope*Udist+intercept+Uene))

        END DO
        DEALLOCATE(Gbin)

      DEALLOCATE(Gbin2)

        CALL RANDOM_NUMBER(rnum)
        SimulTime1=Simultime1+rnum/ratetot
        SimulTime2=SimulTime2-LOG(rnum)/ratetot
        CALL DATE_AND_TIME(datec,timec,zonec,time2)
        CALL TTIME(time1,time2,totaltime)
        subtime1=totaltime
        IF(MOD(jj,1000_int16).EQ.0_int16)THEN
           OPEN(unit=8,file='simul_status.txt',status='UNKNOWN',position='APPEND',action='write')
           WRITE(8,*)jj,naccept,ndepst,subtime1,sfid,nsfatms,nbf,SimulTime1,SimulTime2
          close(8)
        END IF
   END DO   

      DEALLOCATE(xx,xxtemp,ngbharr,ngbharr_br,ngbh1arr_br,fp,f2p,sig,eps,EneArr,MinLocEne)  
! 15   FORMAT(I5,2X,I5,F20.12,2X,F20.12,2X,F20.12,2X,F20.12)
 16   FORMAT(F20.12) 
 17   FORMAT(2X,F20.12,2X,F20.12,2X,F20.12) 
 18   FORMAT(I5,2X,F20.12,2X,F20.12)
      END PROGRAM KMC_DIST_TIMING
