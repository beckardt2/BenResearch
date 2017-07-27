module main

  IMPLICIT NONE
  SAVE

  INTEGER,parameter :: r8=SELECTED_REAL_KIND(16)
  INTEGER,parameter :: int16=SELECTED_INT_KIND(12)
  INTEGER(KIND=int16) :: natms,natms2,frzatms1,mapsiz !Initialized in driver program 
  INTEGER(KIND=int16) ::  intnatms,brnatms,intbrnatms,br1natms,numsteps
  REAL(KIND=r8) :: boxl, rboxl,dxss,dys2,rcut,rcut2,rrcut2,rrcut6,FmaxTol,frzymin,latdist2
  REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: f2p,NewR
  REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: fp,sig,eps,EneArr
  INTEGER(KIND=int16),ALLOCATABLE,DIMENSION(:) :: map,LIST,HEAD,LISTT,HEADT,ngbharr,ngbharr_br,ngbh1arr_br
  INTEGER(KIND=int16),ALLOCATABLE,DIMENSION(:,:) ::COORDNUM 


CONTAINS

  SUBROUTINE conj_grad(x,tol,jmax,Enediff)

    implicit none
     REAL(KIND=r8),INTENT(INOUT) ::x(natms2)
     REAL(KIND=r8),INTENT(IN) :: tol,Enediff
     INTEGER(KIND=int16),INTENT(IN) ::jmax

     INTEGER(KIND=int16) :: j,k,i,imax
     REAL(KIND=r8) :: Uene,deln,beta,delold,delno,delU
     REAL(KIND=r8) :: alp,nfp,nfpold,fno,Uold,Uinit,tol4,xtol,maxd
     REAL(KIND=r8),DIMENSION(natms2) :: r,d,rold,xad,ad
     REAL(KIND=r8),DIMENSION(natms2) :: xold

     x(1:natms)=x(1:natms)-boxl*ANINT(x(1:natms)*rboxl)
     imax=10000_int16
     tol4=tol*tol*tol*tol
     xtol=1e-08_r8

         k=0
         xold=x
         CALL HESSIAN(Uene,x)
         Uinit=Uene

         r=-fp
         d = r
         deln = dot_product(r,r)
         nfp=maxval(abs(fp))
         fno=tol*nfp

         delU=0.0_r8         
           i=0_int16
           DO while(nfp > FmaxTol)
                i=i+1_int16; j=0_int16

              ! DO  while( (j < jmax) .AND. (nfp > fno) )
                DO  while( (j < jmax) .AND. (nfp > FmaxTol) )
                    Uold=Uene
                    alp=-(dot_product(fp,d))/(dot_product(d,matmul(f2p,d)))
                    x=x+alp*d
                    CALL HESSIAN(Uene,x)
                    rold=r
                    nfp=maxval(abs(fp))
                    delold=deln    
                    r=-fp
                    deln=dot_product(r,r)
                    delno=dot_product(r,rold)
                    beta=MAX((deln-delno)/delold,0.0_r8)
                    d=r+beta*d
                    IF( dot_product(r,d) <= 0.0_r8 )d=r
                    delU=Uinit-Uene
                    j=j+1
                END DO
                k=k+1
                if ( (k .EQ. natms2) .OR. dot_product(r,d) <= 0.0_r8) THEN
                    d=r
                    k=0
                END IF

           END DO
           x(1:natms)=x(1:natms)-boxl*ANINT(x(1:natms)*rboxl)

  END SUBROUTINE conj_grad 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE deposit_local_relax(inum2,jmax,xx,tol)
    implicit none

     INTEGER(KIND=int16),INTENT(IN) :: inum2,jmax
     REAL(KIND=r8),INTENT(INOUT) ::xx(natms2)
     REAL(KIND=r8),INTENT(IN) :: tol

!local variables
     INTEGER(KIND=int16) :: inum,imax,i,j,k,ii,n,nx,ny,n2,mm
     REAL(KIND=r8) :: Uene,Uold,deln,delno,nfp,nfpold,fno,beta,alp,delold,tol2,maxd,xtol
     REAL(KIND=r8),DIMENSION(inum2) :: xt,xtold,fsp,r,rold,d,xtad,ad
     REAL(KIND=r8),DIMENSION(inum2,inum2) :: f2sp
     n=natms2; n2=natms; inum=inum2/2
     k=0_int16
     tol2=tol*tol
     xtol=1e-12_r8
!Assign neihboring particles
         DO ii=1,inum
              nx=ngbharr(ii); ny=nx+n2
              xt(ii)=xx(nx); xt(ii+inum)=xx(ny)
         END DO
       
         xtold=xt
         CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xt,xx)
        
         r=-fsp
         d = r
         deln = dot_product(r,r)

         nfp=maxval(abs(fsp))
         fno=tol*nfp
         i=0_int16


           DO while (nfp > FmaxTol) 
                i=i+1_int16; j=0_int16
               DO  while( (j < jmax) .AND. (nfp > FmaxTol) )
                    Uold=Uene
                    alp=-(dot_product(fsp,d))/(dot_product(d,matmul(f2sp,d)))
                    ad=alp*d
                    maxd=maxval(abs(ad))
                    xtad=xt+ad
                    rold=r
                    CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xtad,xx)
                    DO WHILE( ((Uene > Uold) .OR. (maxd > latdist2)) .AND. (maxd > xtol) )
                       alp=0.5_r8*alp
                       ad=alp*d
                       maxd=maxval(abs(ad))
                       xtad=xt+ad
                       CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xtad,xx)
                    END DO
                    xt=xtad

                    r=-fsp
                    nfp=maxval(abs(fsp)) 
                    delold=deln
                    deln=dot_product(r,r)
                    delno=dot_product(r,rold)
                    !beta=MIN((deln-delno)/delold,deln/delold)
                    !beta=MAX(beta,0.0_r8)
                    beta=MAX((deln-delno)/delold,0.0_r8)
                    ! beta=deln/delold
                    d=r+beta*d
                    IF( (Uene > Uold) .OR. (dot_product(r,d) <= 0.0_r8) )d=r
                    j=j+1
                END DO
                k=k+1
                if ( (k .EQ. inum2) .OR. dot_product(r,d) <= 0.0_r8) THEN
                    d=r
                    k=0_int16
                END IF
           END DO

         xt(1:inum)=xt(1:inum)-boxl*ANINT(xt(1:inum)*rboxl)
         DO ii=1,inum
              nx=ngbharr(ii); ny=nx+n2
              xx(nx)=xt(ii); xx(ny)=xt(ii+inum)
         END DO
  END SUBROUTINE deposit_local_relax 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE conj_grad_subset(inum2,jmax,xx,tol)
    implicit none

     INTEGER(KIND=int16),INTENT(IN) :: inum2,jmax
     REAL(KIND=r8),INTENT(INOUT) ::xx(natms2)
     REAL(KIND=r8),INTENT(IN) :: tol

!local variables
     INTEGER(KIND=int16) :: inum,imax,i,j,k,ii,n,nx,ny,n2,mm
     REAL(KIND=r8) :: Uene,Uold,deln,delno,nfp,nfpold,fno,beta,alp,delold,tol2,maxd,xtol
     REAL(KIND=r8),DIMENSION(inum2) :: xt,xtold,fsp,r,rold,d,xtad,ad
     REAL(KIND=r8),DIMENSION(inum2,inum2) :: f2sp
     n=natms2; n2=natms; inum=inum2/2
     k=0_int16; imax=1000_int16; 
     tol2=tol*tol
     xtol=1e-08_r8

!Assign neihboring particles
         DO ii=1,inum
              nx=ngbharr(ii); ny=nx+n2
              xt(ii)=xx(nx); xt(ii+inum)=xx(ny)
         END DO
       
         xtold=xt
         CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xt,xx)
        
         r=-fsp
         d = r
         deln = dot_product(r,r)

         nfp=maxval(abs(fsp))
         fno=tol*nfp
         i=0_int16


           DO while  ( (nfp > FmaxTol) .AND. (i < imax) )
                i=i+1_int16; j=0_int16
               DO  while( (j < jmax) .AND. (nfp > FmaxTol) )
                    Uold=Uene
                    alp=-(dot_product(fsp,d))/(dot_product(d,matmul(f2sp,d)))
                    ad=alp*d
                    maxd=maxval(abs(ad))
                    xtad=xt+ad
                    rold=r
                    CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xtad,xx)
                    DO WHILE( ((Uene > Uold) .OR. (maxd > latdist2)) .AND. (maxd > xtol) )
                       alp=0.5_r8*alp
                       ad=alp*d
                       maxd=maxval(abs(ad))
                       xtad=xt+ad
                       CALL HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xtad,xx)
                    END DO
                    xt=xtad
                    r=-fsp
                    nfp=maxval(abs(fsp)) 
                    delold=deln
                    deln=dot_product(r,r)
                    delno=dot_product(r,rold)
                    !beta=MIN((deln-delno)/delold,deln/delold)
                    !beta=MAX(beta,0.0_r8)
                    beta=MAX((deln-delno)/delold,0.0_r8)
                    ! beta=deln/delold
                    d=r+beta*d
                    IF( (Uene > Uold) .OR. (dot_product(r,d) <= 0.0_r8) )d=r
                    j=j+1

                END DO
                k=k+1
                if ( (k .EQ. inum2) .OR. dot_product(r,d) <= 0.0_r8) THEN
                    d=r
                    k=0_int16
                END IF
           END DO

         xt(1:inum)=xt(1:inum)-boxl*ANINT(xt(1:inum)*rboxl)
         DO ii=1,inum
              nx=ngbharr(ii); ny=nx+n2
              xx(nx)=xt(ii); xx(ny)=xt(ii+inum)
         END DO
  END SUBROUTINE conj_grad_subset 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HESSIAN_SUBSET(inum2,Uene,fsp,f2sp,xt,xx)

    implicit none
     INTEGER(KIND=int16),INTENT(IN) :: inum2
     REAL(KIND=r8),INTENT(INOUT) :: Uene,fsp(inum2),f2sp(inum2,inum2),xt(inum2)
     REAL(KIND=r8),INTENT(IN) :: xx(natms2)

!local variables
    INTEGER(KIND=int16) :: i,j,ix,iy,jx,jy,xj,inum
    REAL(KIND=r8) :: si,ei,xij,xi,yi,xij2,yij,yij2,xyij,rrij,rrij2,rrij3,rrij6
    REAL(KIND=r8) :: xrrij2,yrrij2,sij,sij6,srij,eij,temp,temp1,temp2
    REAL(KIND=r8) :: rijc,scrij,dvdr

    inum=inum2/2
    fsp=0.0_r8
    f2sp=0.0_r8
    Uene=0.0_r8

!Free atoms
   DO ix=1,inum
    iy=ix+inum; i=ngbharr(ix) 
    si=sig(i); ei=eps(i); xi=xt(ix); yi=xt(iy)
    DO jx=ix+1,inum
        jy=jx+inum; xj=ngbharr(jx)

        xij=xi-xt(jx)
        xij=xij-boxl*ANINT(xij*rboxl)
        xij2=xij*xij
        yij=yi-xt(jy)
        yij2=yij*yij
     !   rijc=rcut-SQRT(xij2+yij2)
     ! IF(rijc >= 0.0_r8) THEN

        xyij=xij*yij
        rrij2=1.0_r8/(xij2 + yij2)
        rrij=sqrt(rrij2)
        rrij3=rrij2*rrij
        rrij6=rrij3*rrij3
        
        xrrij2=xij2*rrij2  
        yrrij2=yij2*rrij2
        sij=0.5_r8*(si+sig(xj))   
     !   sij6=sij**6
     !   srij=sij6 * rrij6
     !   scrij = sij6*rrcut6

     !   eij=4.0_r8*sqrt(ei*eps(xj))
     !   dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
     !   Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc

     !   !Derivatives of the potential
        
     !   temp=6.0_r8*eij*rrij2*srij*(2.0_r8*srij-1.0_r8)+dvdr
 
        srij=sij**6 * rrij6

        eij=4.0_r8*srij*sqrt(ei*eps(xj))
        Uene=Uene + eij*(srij-1.0_r8)
        !Derivatives of the potential

        temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)

        fsp(ix)=fsp(ix)-temp*xij
        fsp(iy)=fsp(iy)-temp*yij

        fsp(jx)=fsp(jx)+temp*xij
        fsp(jy)=fsp(jy)+temp*yij

        !Hessian matrix
        temp1=6.0_r8*eij*rrij*(1.0_r8-2.0_r8*srij)
        temp2=6.0_r8*eij*rrij2*(26.0_r8*srij-7.0_r8)
        
        !Hessian matrix
     !   temp1=6.0_r8*eij*rrij*srij*(1.0_r8-2.0_r8*srij)
     !   temp2=6.0_r8*eij*rrij2*srij*(26.0_r8*srij-7.0_r8)
        
        f2sp(ix,ix)=f2sp(ix,ix)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
        f2sp(ix,iy)=f2sp(ix,iy)-temp1*xyij*rrij3+temp2*xyij*rrij2
        
        f2sp(ix,jx)=temp1*rrij*(xrrij2-1.0_r8)-temp2*xrrij2
        f2sp(ix,jy)=temp1*xyij*rrij3-temp2*xyij*rrij2
        f2sp(jx,ix)=f2sp(ix,jx)
        f2sp(jy,ix)=f2sp(ix,jy)
        
        f2sp(iy,iy)=f2sp(iy,iy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2
        
        f2sp(iy,jy)=temp1*rrij*(yrrij2-1.0_r8)-temp2*yrrij2
        f2sp(iy,jx)=temp1*xyij*rrij3-temp2*xyij*rrij2  
        f2sp(jy,iy)=f2sp(iy,jy)
        f2sp(jx,iy)=f2sp(iy,jx)
        
        f2sp(jx,jx)=f2sp(jx,jx)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
        f2sp(jx,jy)=f2sp(jx,jy)-temp1*xyij*rrij3+temp2*xyij*rrij2
        f2sp(jy,jy)=f2sp(jy,jy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2
        
        f2sp(jy,jx)=f2sp(jx,jy)
     !  END IF
    END DO            
    f2sp(iy,ix)=f2sp(ix,iy)
    
    DO j=1,brnatms
            jx=ngbharr_br(j)
            jy=jx+natms
            xij=xi-xx(jx)
            xij=xij-boxl*ANINT(xij*rboxl)
            xij2=xij*xij
            yij=yi-xx(jy)
            yij2=yij*yij
      !      rijc=rcut-SQRT(xij2+yij2)
      !  IF(rijc >= 0.0_r8) THEN

            xyij=xij*yij
            rrij2=1.0_r8/(xij*xij + yij*yij)
            xrrij2=xij2*rrij2
            yrrij2=yij2*rrij2

            rrij=sqrt(rrij2)
            rrij3=rrij2*rrij
            rrij6=rrij3*rrij3

            sij=0.5_r8*(si+sig(jx))   
            srij=sij**6 *rrij6

            eij=4.0_r8*srij*sqrt(ei*eps(jx))
            Uene=Uene + eij*(srij-1.0_r8)
            !Derivatives of the potential

            temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)

            !Forces and Hessian matrix
            fsp(ix)=fsp(ix)-temp*xij
            fsp(iy)=fsp(iy)-temp*yij

            !Hessian matrix
            temp1=6.0_r8*eij*rrij*(1.0_r8-2.0_r8*srij)
            temp2=6.0_r8*eij*rrij2*(26.0_r8*srij-7.0_r8)


       !     sij6=sij**6
       !     srij=sij6 *rrij6

       !     scrij = sij6*rrcut6

       ! !Shifted force-potential

       !     eij=4.0_r8*sqrt(ei*eps(xj))
       !     dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
       !     Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc

       !     temp=6.0_r8*eij*rrij2*srij*(2.0_r8*srij-1.0_r8)+dvdr


        !Hessian matrix
       !     temp1=6.0_r8*eij*rrij*srij*(1.0_r8-2.0_r8*srij)
       !     temp2=6.0_r8*eij*rrij2*srij*(26.0_r8*srij-7.0_r8)

            f2sp(ix,ix)=f2sp(ix,ix)+0.5_r8*(temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2)
            f2sp(ix,iy)=f2sp(ix,iy)-0.5_r8*xyij*(temp1*rrij3-temp2*rrij2)
            f2sp(iy,iy)=f2sp(iy,iy)+0.5_r8*(temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2)
            f2sp(iy,ix)=f2sp(ix,iy)

       !   END IF
    END DO

   END DO 

 RETURN           
 END SUBROUTINE HESSIAN_SUBSET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE SEARCH_NGBHARR(jx,inum,xj,test)

   implicit none
     INTEGER(KIND=int16),INTENT(IN) :: jx,inum
     INTEGER(KIND=int16),INTENT(INOUT) :: xj
     LOGICAL,INTENT(INOUT) :: test

!local variables
     INTEGER(KIND=int16) :: i

     test=.FALSE.
     xj=jx
     DO i=1,inum
        IF(jx==ngbharr(i)) THEN
           xj=i;
           TEST=.TRUE.
           EXIT
        END IF
     END DO 
 
 END SUBROUTINE SEARCH_NGBHARR 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE STEEP_DESC(ix,x,tol,nsteps)

    implicit none
     INTEGER(KIND=int16),INTENT(IN) :: ix
     INTEGER(KIND=int16),INTENT(OUT) :: nsteps
     REAL(KIND=r8),INTENT(INOUT) ::x(natms2)
     REAL(KIND=r8),INTENT(IN) :: tol

!local variables
     INTEGER(KIND=int16) :: iy,i,imax
     REAL(KIND=r8), PARAMETER :: eps0=0.01_r8
     REAL(KIND=r8) :: xxold(2),uold,uene,eps,ftol,fno,xx(2),fpo(2),nfpo

     ftol=tol*FmaxTol
     imax=10000_int16
     eps=eps0
     iy=ix+natms
     xx(1)=x(ix); xx(2)=x(iy)
     CALL POT_FORCE_ON_ONE(ix,x,uene,fpo)
     nfpo=maxval(abs(fpo))
     fno=tol*nfpo
     i=0_int16
!     DO while ( (maxval(abs(fpo)) > fno) .AND. (i .LE. imax))
     DO while ( (nfpo > fno) .AND. (nfpo > fTol) )
           i=i+1
           uold=uene
           xxold=xx;
           xx=xxold-eps*fpo
           x(ix)=xx(1); x(iy)=xx(2)
           CALL POT_FORCE_ON_ONE(ix,x,uene,fpo)
           IF(uene > uold) THEN
             xx=xxold; eps=0.5_r8*eps 
             CYCLE
           END IF
           eps=eps0
           nfpo=maxval(abs(fpo))
           x(ix)=x(ix)-boxl*ANINT(x(ix)*rboxl)
     END DO
     nsteps=i
 
 END SUBROUTINE STEEP_DESC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE COMPUTE_REF_DIST(midrow,natmsw,x,dxs,dHref,dDref)
     
    implicit none
     INTEGER(KIND=int16),INTENT(IN) :: midrow,natmsw
     REAL(KIND=r8),INTENT(IN) :: x(natms2),dxs
     REAL(KIND=r8),INTENT(OUT) :: dHref,dDref

!local variables
     INTEGER(KIND=int16) :: ix,iy,jx,jy,j,numH,numD,Abin(1,2)
     REAL(KIND=r8) :: xi,yi,xij,yij,dxd,dxd_br

     dHref=0.0_r8; dDref=0.0_r8; numH=0_int16; numD=0_int16
     DO ix=midrow,midrow+natmsw
        Abin(1,1)=ix; Abin(1,2)=0_int16
        dxd=dxs; dxd_br=dxd+1.0_r8*dxs
        CALL COORD_NUM_SIMPLE(1_int16,Abin,dxd,dxd_br,x)

        iy=ix+natms; xi=x(ix); yi=x(iy)

        DO j=2,intnatms
          jx=ngbharr(j)
          jy=jx+natms

          xij=xi-x(jx); yij=ABS(yi-x(jy))
          xij=xij-boxl*ANINT(xij*rboxl)
          IF(yij > 1.0_r8) THEN
             dDref = dDref + SQRT(xij*xij + yij*yij)
             numD  = numD  + 1_int16
          ELSE
             dHref = dHref + SQRT(xij*xij + yij*yij) 
             numH  = numH  + 1_int16
          END IF
        END DO

     END DO
     dDref = dDref/REAL(numD,kind=r8)
     dHref = dHref/REAL(numH,kind=r8)

  END SUBROUTINE COMPUTE_REF_DIST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE COMPUTE_STRAIN(Mx,My,dHref,dDref,dxs,x,strain,strainxy)

    implicit none
     INTEGER(KIND=int16),INTENT(IN) :: Mx,My
     REAL(KIND=r8),INTENT(IN) :: dxs,dHref,dDref
     REAL(KIND=r8),INTENT(INOUT) :: x(natms2)
     REAL(KIND=r8),INTENT(OUT) :: strain(natms),strainxy(natms,2)

!local variables
     INTEGER(KIND=int16) :: ix,iy,jx,jy,j,Abin(1,2)
     REAL(KIND=r8) :: xi,yi,xij,yij,dxd,dxd_br,disp

     strain=0.0_r8; strainxy=0.0_r8
     CALL LINKS(Mx,My,x)
     DO ix=1,natms
        Abin(1,1)=ix; Abin(1,2)=0_int16
        dxd=dxs; dxd_br=dxd+1.0_r8*dxs
        CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,x)

        iy=ix+natms; xi=x(ix); yi=x(iy)

        DO j=2,intnatms
          jx=ngbharr(j)
          jy=jx+natms

          xij=xi-x(jx); yij=ABS(yi-x(jy))
          xij=xij-boxl*ANINT(xij*rboxl)
          IF(yij > 1.0_r8) THEN
             disp = SQRT(xij*xij + yij*yij) - dDref
             strain(ix) = strain(ix) + disp*disp
             strainxy(ix,2)=strainxy(ix,2)+disp*disp
          ELSE
             disp = SQRT(xij*xij + yij*yij) - dHref
             strain(ix) = strain(ix) + disp*disp
             strainxy(ix,1)=strainxy(ix,1)+disp*disp
          END IF
        END DO

     END DO
     Strain=SQRT(Strain) 

  END SUBROUTINE COMPUTE_STRAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE POT_FORCE_ON_ONE(ix,x,uene,fpo)

    implicit none
     INTEGER(KIND=int16),INTENT(IN) :: ix
     REAL(KIND=r8),INTENT(INOUT) ::x(natms2),uene,fpo(2)

!local variables
     INTEGER(KIND=int16) :: j,iy,jx,jy
     REAL(KIND=r8) :: si,ei,sij,sij6,srij,eij,xi,yi,xij,yij,rrij2,rrij6,temp
     REAL(KIND=r8) :: dvdr,scrij,rijc
     uene=0.0_r8    
     fpo=0.0_r8

     iy=ix+natms; xi=x(ix); yi=x(iy); si=sig(ix); ei=eps(ix)

    DO j=2,intnatms
        jx=ngbharr(j)
        jy=jx+natms

        xij=xi-x(jx); yij=yi-x(jy)
        xij=xij-boxl*ANINT(xij*rboxl)
       ! rijc=rcut-SQRT(xij*xij+yij*yij)
      !IF(rijc >= 0.0_r8) THEN

        rrij2=1.0_r8/(xij*xij+yij*yij)
        rrij6=rrij2*rrij2*rrij2

        sij=0.5_r8*(si+sig(jx))
        srij=sij**6 * rrij6

        eij=4.0_r8*srij*SQRT(ei*eps(jx))
        uene=uene + eij*(srij-1.0_r8)

        temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)
        fpo(1)=fpo(1)-temp*xij
        fpo(2)=fpo(2)-temp*yij

     !   sij6=sij**6
     !   srij=sij6 * rrij6
     !   scrij=sij6*rrcut6

     !   eij=4.0_r8*SQRT(ei*eps(jx))
     !   dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
     !   Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc

     !   temp=6.0_r8*eij*rrij2*srij*(2.0_r8*srij-1.0_r8)+dvdr
      !END IF
    END DO

    DO j=1,brnatms
        jx=ngbharr_br(j)
        jy=jx+natms

        xij=xi-x(jx); yij=yi-x(jy)
        xij=xij-boxl*ANINT(xij*rboxl)
      !  rijc=rcut-SQRT(xij*xij+yij*yij)
      !IF(rijc >= 0.0_r8) THEN

        rrij2=1.0_r8/(xij*xij+yij*yij)
        rrij6=rrij2*rrij2*rrij2

        sij=0.5_r8*(si+sig(jx))
        srij=sij**6 * rrij6

        eij=4.0_r8*srij*SQRT(ei*eps(jx))
        uene=uene + eij*(srij-1.0_r8)

        temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)
        fpo(1)=fpo(1)-temp*xij
        fpo(2)=fpo(2)-temp*yij

     !   sij6=sij**6
     !   srij=sij6 * rrij6
     !   scrij=sij6*rrcut6

     !   eij=4.0_r8*SQRT(ei*eps(jx))
     !   dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
     !   Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc

     !   temp=6.0_r8*eij*rrij2*srij*(2.0_r8*srij-1.0_r8)+dvdr
     ! END IF

    END DO


 END SUBROUTINE POT_FORCE_ON_ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE BOUNDARY_FORCE(xx,nbf,intst)
    implicit none
     REAL(KIND=r8),INTENT(OUT) :: nbf 
     REAL(KIND=r8),INTENT(IN) :: xx(natms2)
     LOGICAL :: intst
!local variables
    INTEGER(KIND=int16) :: i,i1,j,ix,iy,jx,jy
    REAL(KIND=r8) :: si,ei,xij,xi,yi,xij2,yij,yij2,rrij2,rrij6
    REAL(KIND=r8) :: sij,sij6,srij,eij,temp,fsp(2_int16*br1natms)
    REAL(KIND=r8) :: dvdr,scrij,rijc,Uene

    intst=.FALSE.
    fsp=0.0_r8; Uene=0.0_r8
 DO i=1,br1natms
    ix=ngbh1arr_br(i); iy=ix+natms
    xi=xx(ix); yi=xx(iy)
    si=sig(ix); ei=eps(ix)
    i1=i+br1natms
 
    DO j=1,brnatms

      jx=ngbharr_br(j); 

      IF(ix .NE. jx) THEN
        jy=jx+natms


        xij=xi-xx(jx)
        xij=xij-boxl*ANINT(xij*rboxl)
        xij2=xij*xij
        yij=yi-xx(jy)
        yij2=yij*yij

        rrij2=1.0_r8/(xij2+yij2)
        rrij6=rrij2*rrij2*rrij2

        sij=0.5_r8*(si+sig(jx))
        srij=sij**6 * rrij6

        eij=4.0_r8*srij*sqrt(ei*eps(jx))
        Uene=Uene + eij*(srij-1.0_r8)
        !Derivatives of the potential

        temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)

        fsp(i)=fsp(i)-temp*xij
        fsp(i1)=fsp(i1)-temp*yij

      END IF
    END DO

    DO j=1,intnatms
        jx=ngbharr(j); jy=jx+natms

        xij=xi-xx(jx)
        xij=xij-boxl*ANINT(xij*rboxl)
        xij2=xij*xij
        yij=yi-xx(jy)
        yij2=yij*yij

        rrij2=1.0_r8/(xij2+yij2)
        rrij6=rrij2*rrij2*rrij2

        sij=0.5_r8*(si+sig(jx))
        srij=sij**6 * rrij6

        eij=4.0_r8*srij*sqrt(ei*eps(jx))
        Uene=Uene + eij*(srij-1.0_r8)
        !Derivatives of the potential

        temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)

        fsp(i)=fsp(i)-temp*xij
        fsp(i1)=fsp(i1)-temp*yij


        IF(xx(jy) < frzymin)intst=.TRUE.
          
    END DO

  END DO

  nbf=MAXVAL(abs(fsp))

 END SUBROUTINE BOUNDARY_FORCE
!!!
 SUBROUTINE HESSIAN_CUTOFF(Uene,x)

   IMPLICIT NONE

    REAL(KIND=r8),INTENT(INOUT) :: Uene
    REAL(KIND=r8),INTENT(IN) :: x(natms2)

    INTEGER(KIND=int16) :: ix,iy,j,jx,jy,m1,m2,len,len2,ll,Abin(1,2)
    REAL(KIND=r8) :: xix,yix,six,eix,xij,xij2,yij,yij2,xyij,rrij,rrij2,rrij3,rrij6
    REAL(KIND=r8) :: xrrij2,yrrij2,sij,srij,eij,temp,temp1,temp2,dxs,dxd,dxd_br
    REAL(KIND=r8) :: sigt(natms),epst(natms)

    sigt=sig; epst=eps
    dxs=2.318_r8

    fp=0.0_r8
    f2p=0.0_r8
    Uene=0.0_r8

    DO ix=1,natms
       iy=ix+natms; xix=x(ix); yix=x(iy); six=sig(ix); eix=eps(ix);
       dxd=5.0_r8*dxs; dxd_br=dxd+dxs;
       Abin(1,1)=ix; Abin(1,2)=0_int16
       CALL COORD_NUM_SIMPLE(1_int16,Abin,dxd,dxd_br,x)

       DO j=2,intnatms
          jx=ngbharr(j)
          jy=jx+natms

          xij=xix-x(jx); xij=xij-boxl*ANINT(xij*rboxl)
          yij=yix-x(jy)
          xij2=xij*xij
          yij2=yij*yij
          xyij=xij*yij
          rrij2=1.0_r8/(xij*xij + yij*yij)
          rrij=sqrt(rrij2)
          rrij3=rrij2*rrij
          rrij6=rrij3*rrij3

          xrrij2=xij2*rrij2
          yrrij2=yij2*rrij2
          sij=0.5_r8*(six+sigt(jx))
          srij=sij**6*rrij6

          eij=4.0_r8*srij*sqrt(eix*epst(jx))
          Uene=Uene + eij*(srij-1.0_r8)
        !Derivatives of the potential
        
          temp=6.0_r8*eij*rrij2*(2.0_r8*srij-1.0_r8)
        
          fp(ix)=fp(ix)-temp*xij
          fp(iy)=fp(iy)-temp*yij
          
          fp(jx)=fp(jx)+temp*xij
          fp(jy)=fp(jy)+temp*yij   
        !Hessian matrix
          temp1=6.0_r8*eij*rrij*(1.0_r8-2.0_r8*srij)
          temp2=6.0_r8*eij*rrij2*(26.0_r8*srij-7.0_r8)
        
          f2p(ix,ix)=f2p(ix,ix)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
          f2p(ix,iy)=f2p(ix,iy)-temp1*xyij*rrij3+temp2*xyij*rrij2
        
          f2p(ix,jx)=temp1*rrij*(xrrij2-1.0_r8)-temp2*xrrij2
          f2p(ix,jy)=temp1*xyij*rrij3-temp2*xyij*rrij2
          f2p(jx,ix)=f2p(ix,jx)
          f2p(jy,ix)=f2p(ix,jy)
        
          f2p(iy,iy)=f2p(iy,iy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2
        
          f2p(iy,jy)=temp1*rrij*(yrrij2-1.0_r8)-temp2*yrrij2
          f2p(iy,jx)=temp1*xyij*rrij3-temp2*xyij*rrij2  
          f2p(jy,iy)=f2p(iy,jy)
          f2p(jx,iy)=f2p(iy,jx)
        
          f2p(jx,jx)=f2p(jx,jx)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
          f2p(jx,jy)=f2p(jx,jy)-temp1*xyij*rrij3+temp2*xyij*rrij2
          f2p(jy,jy)=f2p(jy,jy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2

          f2p(jy,jx)=f2p(jx,jy)
 
       END DO  
       sigt(ix)=0.0_r8; epst(ix)=0.0_r8
       f2p(iy,ix)=f2p(ix,iy)    
   END DO
   

END SUBROUTINE HESSIAN_CUTOFF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 SUBROUTINE HESSIAN(Uene,x)

   IMPLICIT NONE
    
    REAL(KIND=r8),INTENT(INOUT) :: Uene
    REAL(KIND=r8),INTENT(IN) :: x(natms2)

    INTEGER(KIND=int16) :: ix,iy,jx,jy,m1,m2,len,len2,ll
    REAL(KIND=r8) :: xix,yix,six,eix,xij,xij2,yij,yij2,xyij,rrij,rrij2,rrij3,rrij6
    REAL(KIND=r8) :: xrrij2,yrrij2,sij,sij6,srij,eij,temp,temp1,temp2
    REAL(KIND=r8) :: rijc,scrij,dvdr

    fp=0.0_r8
    f2p=0.0_r8
    Uene=0.0_r8
    
    DO ix=1,natms-1
       iy=ix+natms; xix=x(ix); yix=x(iy); six=sig(ix); eix=eps(ix);
       DO jx=ix+1,natms
          jy=jx+natms

          xij=xix-x(jx); xij=xij-boxl*ANINT(xij*rboxl)
          yij=yix-x(jy)
          xij2=xij*xij
          yij2=yij*yij
          rijc=rcut-SQRT(xij2+yij2)
      IF(rijc >= 0.0_r8) THEN
          xyij=xij*yij
          rrij2=1.0_r8/(xij*xij + yij*yij)
          rrij=sqrt(rrij2)
          rrij3=rrij2*rrij
          rrij6=rrij3*rrij3
        
          xrrij2=xij2*rrij2  
          yrrij2=yij2*rrij2
          sij=0.5_r8*(six+sig(jx))   
          sij6=sij**6
          srij=sij6*rrij6
          scrij = sij6*rrcut6
  
          eij=4.0_r8*sqrt(eix*eps(jx))
          dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
          Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc
        !Derivatives of the potential
        
          temp=6.0_r8*eij*rrij2*srij*(2.0_r8*srij-1.0_r8)+dvdr
        
          fp(ix)=fp(ix)-temp*xij
          fp(iy)=fp(iy)-temp*yij
          
          fp(jx)=fp(jx)+temp*xij
          fp(jy)=fp(jy)+temp*yij   

        !Hessian matrix
          temp1=6.0_r8*eij*rrij*srij*(1.0_r8-2.0_r8*srij)
          temp2=6.0_r8*eij*rrij2*srij*(26.0_r8*srij-7.0_r8)
        
          f2p(ix,ix)=f2p(ix,ix)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
          f2p(ix,iy)=f2p(ix,iy)-temp1*xyij*rrij3+temp2*xyij*rrij2
        
          f2p(ix,jx)=temp1*rrij*(xrrij2-1.0_r8)-temp2*xrrij2
          f2p(ix,jy)=temp1*xyij*rrij3-temp2*xyij*rrij2
          f2p(jx,ix)=f2p(ix,jx)
          f2p(jy,ix)=f2p(ix,jy)
       
          f2p(iy,iy)=f2p(iy,iy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2
        
          f2p(iy,jy)=temp1*rrij*(yrrij2-1.0_r8)-temp2*yrrij2
          f2p(iy,jx)=temp1*xyij*rrij3-temp2*xyij*rrij2  
          f2p(jy,iy)=f2p(iy,jy)
          f2p(jx,iy)=f2p(iy,jx)
        
          f2p(jx,jx)=f2p(jx,jx)+temp1*rrij*(1.0_r8-xrrij2)+temp2*xrrij2
          f2p(jx,jy)=f2p(jx,jy)-temp1*xyij*rrij3+temp2*xyij*rrij2
          f2p(jy,jy)=f2p(jy,jy)+temp1*rrij*(1.0_r8-yrrij2)+temp2*yrrij2

          f2p(jy,jx)=f2p(jx,jy)
        END IF
     END DO            
       f2p(iy,ix)=f2p(ix,iy)    
   END DO
   
END SUBROUTINE HESSIAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SETUP_MAP(mx,my)
     
      IMPLICIT NONE

      INTEGER(KIND=int16),INTENT(IN) :: mx,my

! local variables
      INTEGER(KIND=int16),DIMENSION(3,26) :: ngbh_list
      INTEGER(KIND=int16) ::  ix,iy,iz,imap,index


      ngbh_list=0_int16
! create ngbh_list with 26 direction vectors (in 3D)
!or  create ngbh_list with 8 direction vectors (in 2D)

      index=0_int16
      DO ix=-1,1
         DO iy=-1,1
!            DO iz=-1,1
!               IF ((ix*ix+iy*iy+iz*iz) .NE. 0) THEN
               IF ((ix*ix+iy*iy) .NE. 0) THEN
                  index=index+1
                  ngbh_list(1,index)=ix
                  ngbh_list(2,index)=iy
!                  ngbh_list(3,index)=iz
               END IF
!            END DO
         END DO
      END DO 


!    find the 26 neighbors of each cell (in 3D)
!    find the 8 neighbors of each cell (in 2D)

!      DO iz = 1, m
         DO iy = 1, my
            DO ix = 1, mx

!               imap = (cell(ix,iy,iz,m)-1)*26
!               DO index=1,26
!                  map(imap+index)=cell(ix+ngbh_list(1,index),&
!                                      iy+ngbh_list(2,index),&
!                                      iz+ngbh_list(3,index),m)
!               END DO

               imap = (cell(ix,iy,mx,my)-1)*8
               DO index=1,8
                  map(imap+index)=cell(ix+ngbh_list(1,index),&
                                      iy+ngbh_list(2,index),mx,my)
               END DO

            END DO
         END DO
!      END DO
      
      END SUBROUTINE SETUP_MAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER(KIND=int16) FUNCTION CELL(ix,iy,mx,my)
!      INTEGER(KIND=int16) FUNCTION CELL(ix,iy,iz,m)
      IMPLICIT NONE
                                                                
      INTEGER(KIND=int16),INTENT(IN) :: ix,iy,mx,my !,iz
                                                               
!      cell = 1 + MOD(ix-1+m,m) &
!                + MOD(iy-1+m,m) * m &
!                + MOD(iz-1+m,m) * m*m 
                                                  
      cell = 1 + MOD(ix-1+mx,mx) &
                + MOD(iy-1+my,my) * mx 
             
      RETURN
      END FUNCTION CELL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(KIND=r8) FUNCTION BOUNDS(misfit)
        IMPLICIT NONE

        REAL(KIND=r8),INTENT(IN) :: misfit
        REAL(KIND=r8),PARAMETER  :: aa=0.008815832398868_r8
        REAL(KIND=r8),PARAMETER  :: bb=0.108611252923514_r8
        REAL(KIND=r8),PARAMETER  :: cc=1.626164470437413_r8
        REAL(KIND=r8) :: msft 

          msft=100.0_r8*misfit
          bounds=aa*msft*msft + bb*msft + cc
      RETURN
      END FUNCTION BOUNDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(KIND=r8) FUNCTION BOUNDS_BIEHL(misfit)
        IMPLICIT NONE

        REAL(KIND=r8),INTENT(IN) :: misfit
        REAL(KIND=r8),PARAMETER  :: aa=-0.023898129637131_r8
        REAL(KIND=r8),PARAMETER  :: bb=0.323364990861612_r8
        REAL(KIND=r8),PARAMETER  :: cc=1.190093851255458_r8
       REAL(KIND=r8) :: msft 
 
          msft=100.0_r8*misfit
          bounds_biehl=aa*msft*msft + bb*msft + cc
      RETURN
      END FUNCTION BOUNDS_BIEHL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LINKS_ONE(ix,Mx,My,xxtemp,xx)

!This routine adusts the LINKED list when the position of just one atom
!changes
   
      IMPLICIT NONE
        INTEGER(KIND=int16),INTENT(IN) ::  ix,Mx,My
        REAL(KIND=r8),INTENT(IN) :: xxtemp(natms2),xx(natms2)

!local variables
        INTEGER(KIND=int16) :: iy,jx,jjx,iix,frz1,icello,icelln
        REAL(KIND=r8) :: rboxly,cellix,celliy!,celliz

        iy=ix+natms
        rboxly=1.0_r8/(MAXVAL(xxtemp(natms+1:natms2))+dys2)

        cellix=REAL(Mx,KIND=r8); celliy=REAL(My,KIND=r8)

!Find atom's old cell
        icello = 1_int16 + INT((xxtemp(ix)*rboxl + 0.5_r8 ) * cellix ) & 
                         + INT( xxtemp(iy)*rboxly           * celliy ) * Mx 

        rboxly=1.0_r8/(MAXVAL(XX(natms+1:natms2))+dys2)
!Find atom's new cell
        icelln = 1_int16 + INT((xx(ix)*rboxl + 0.5_r8 ) * cellix ) & 
                         + INT( xx(iy)*rboxly           * celliy ) * Mx 

        IF(icelln .NE. icello) THEN ! Atom has new cell, adjust list accordingly
          
           IF(HEAD(icello)==ix) THEN !If atom is stored in head then
 
              HEAD(icello)=list(ix)  !store the next atom in the list in head
           ELSE   !Atom is not in head, it's in LIST, so adjust list by removing
                  !atom 

              jx=head(icello); jx=list(jx); jjx=jx
              IF(jx == ix) THEN  !If atom is the first atom in list, then
                                 !just point head to the second atom in list
                 HEAD(icello)=jx
              ELSE  !Atom is the not the first in the list, so we have to
                    !remove atom from this part of the list
                 DO WHILE (jx > 0_int16)

                    jx=list(jx) !Find the next atom in the list
                    IF( jx == ix) THEN    !If next atom is what we are looking for then
                      list(jjx) = list(jx) !point previous atom to the atom after iix 
                      EXIT                 ! and exit do loop
                    END IF
                    jjx=jx  !If not, then set jjx to present atom and repeat loop
                     
                 END DO

              END IF

            END IF

         !Now insert atom into its new list; We will insert it into the last spot which stores 0

            jx=head(icelln); jjx=jx; jx=list(jx)                    
            DO WHILE (jx > 0_int16)
               jjx=jx
               jx=list(jx)
            END DO
            list(jjx)=ix
            list(ix)=0_int16

        END IF
        END SUBROUTINE LINKS_ONE
!!!!!!!!!!!!!!!!!!
      SUBROUTINE LINKS(Mx,My,XX)
    
      IMPLICIT NONE
        INTEGER(KIND=int16),INTENT(IN) ::  Mx,My
        REAL(KIND=r8),INTENT(INOUT) :: XX(natms2)

!local variables
        INTEGER(KIND=int16) :: i,err,icell,xcell,ycell,zcell
        REAL(KIND=r8) :: RX(natms),RY(natms),rboxly,cellix,celliy,miny,maxy,aminy!,celliz

        XX(1:natms)=XX(1:natms)-boxl*ANINT(XX(1:natms)*rboxl)
        RX=XX(1:natms); RX=RX*rboxl
        RY=XX(natms+1:natms2); miny=MINVAL(RY); maxy=MAXVAL(RY); aminy=abs(miny)
        rboxly=1.0_r8/(maxy-miny+dys2)
        RY=RY+aminy

        cellix=REAL(Mx,KIND=r8); celliy=REAL(My,KIND=r8)*rboxly

        HEAD=0_int16; LIST=0_int16;

        DO I = 1, natms
           
           ICELL = 1_int16 + INT ( ( RX(I) + 0.5_r8 ) * cellix ) & ! -L/2 <= x <= L/2
                           + INT (   RY(I)            * celliy ) * Mx !&
!                           + INT ( ( RZ(I) + 0.5_r8 ) * celliz ) * M * M
          
           xcell = INT ( ( RX(I) + 0.5_r8 ) * cellix ) 
           ycell = INT ( RY(I) * celliy ) * Mx
!           ICELL = 1_int16 + xcell + ycell 
            IF((ICELL > MX*MY) .OR. (ICELL <= 0)) THEN
               WRITE(4,*)i,xcell,ycell,icell,xx(i),ry(i)
               WRITE(4,*) boxl,rboxl,rboxly
            END IF
           LIST(I)     = HEAD(ICELL)
           HEAD(ICELL) = I
        END DO

        
        END SUBROUTINE LINKS
!!!!!!!!!!!!!!!!!!
      SUBROUTINE COORD_NUM(GBLEN,Gbin,Mx,My,dxd,dxd_br,XX)

      IMPLICIT NONE
        INTEGER(KIND=int16),INTENT(IN) :: GBLEN,Mx,My
        INTEGER(KIND=int16),INTENT(INOUT) :: Gbin(GBLEN,2)
        REAL(KIND=r8),INTENT(IN) :: dxd,dxd_br,XX(natms2)

!local variables
       INTEGER(KIND=int16) :: i,ix,iy,jx,jy,coord,icell,jcell0,jcell,nabor
       INTEGER(KIND=int16) :: xcell,ycell,zcell
       REAL(KIND=r8) :: rboxly,cellix,celliy,maxy,miny,aminy
       REAL(KIND=r8) :: dxsd,dxsd2,dxd_brn,dxd_brn2,xix,yix,yjx,xij,yij,rij2
       REAL(KIND=r8) :: dxbd,dxbd2
       frzymin=0.0_r8

       dxsd=dxd; dxsd2=dxsd*dxsd
       dxd_brn=dxd_br + 0.2_r8; dxd_brn2=dxd_brn*dxd_brn
       dxbd=dxd+dxss+0.2; dxbd2=dxbd*dxbd
       maxy=MAXVAL(XX(natms+1:natms2)); miny=MINVAL(XX(natms+1:natms2)); aminy=abs(miny)
       rboxly=1.0_r8/(maxy-miny + dys2)
       cellix=REAL(Mx,KIND=r8); celliy=REAL(My,KIND=r8)*rboxly
       
       DO i=1,GBLEN

          coord=0_int16; brnatms=0_int16; br1natms=0_int16
          ix=Gbin(i,1); iy=ix+natms; ngbharr(1)=ix; 
          IF((ix < 1) .OR. (ix > natms))WRITE(3,*)i,ix
          xix=XX(ix); yix=XX(iy)
          IF(ix > 0) THEN
           ICELL = 1 + INT ( ( XX(ix)*rboxl + 0.5_r8 ) * cellix ) &
                     + INT (    (XX(iy)+aminy)       * celliy ) * Mx !&
!                    + INT ( ( RZ(iz) + 0.5_r8 ) * celliz ) * M * M

!           xcell = INT ( ( XX(ix)*rboxl + 0.5_r8 ) * cellix )
!           ycell = MIN(M*(M-1_int16), INT ( XX(iy) * celliy ) * Mx)
!           ICELL = 1_int16 + xcell + ycell
 
           jx=HEAD(ICELL); jy=jx+natms

!Loop over particles in the same cell
           DO WHILE ( jx > 0 )  

             yjx=XX(jy)
             xij=xix-XX(jx); xij=xij-boxl*ANINT(xij*rboxl)
             yij=yix-yjx
             rij2=xij*xij+yij*yij
!             IF( (rij2 > 0.0_r8) .AND. (rij2 <= dxsd2)) THEN
             IF( (jx .NE. ix) .AND. (rij2 <= dxsd2)) THEN
                 coord=coord+1_int16
                 ngbharr(coord+1)=jx
             ELSEIF( (rij2 > dxsd2) .AND. (rij2 <= dxd_brn2)) THEN
                 IF(rij2 <= dxbd2) THEN
                   br1natms=br1natms+1_int16
                   ngbh1arr_br(br1natms)=jx
                 END IF
                 brnatms=brnatms+1_int16
                 ngbharr_br(brnatms)=jx
                 IF(yjx <= frzymin)frzymin=yjx
             END IF
             jx=LIST(jx); jy=jx+natms

           END DO

!Loop over neighboring 26 cells (in 3D)
!           jcell0=26_int16*(icell-1_int16)
!           DO nabor=1,26
             

!Loop over neighboring 8 cells (in 2D)
           jcell0=8_int16*(icell-1_int16)
           DO nabor=1,8
              jcell=map(jcell0+nabor)
          !Loop over all molcules in neighboring cells
              jx=HEAD(JCELL); jy=jx+natms

              DO WHILE ( jx > 0 ) 

                yjx=XX(jy)
                xij=xix-XX(jx); xij=xij-boxl*ANINT(xij*rboxl)
                yij=yix-yjx
                rij2=xij*xij+yij*yij
                IF( (rij2 > 0.0_r8) .AND. (rij2 <= dxsd2)) THEN
                    coord=coord+1_int16
                    ngbharr(coord+1)=jx
                ELSEIF( (rij2 > dxsd2) .AND. (rij2 <= dxd_brn2)) THEN
                    IF(rij2 <= dxbd2) THEN
                       br1natms=br1natms+1_int16
                       ngbh1arr_br(br1natms)=jx
                    END IF
                    brnatms=brnatms+1_int16
                    ngbharr_br(brnatms)=jx
                    IF(yjx <= frzymin)frzymin=yjx
                END IF

                jx=LIST(jx); jy=jx+natms

             END DO

           END DO
           Gbin(i,2)=coord
          END IF

        END DO
        intnatms=coord+1
       
        END SUBROUTINE COORD_NUM
!!!!!!!!!!!!!!!!!!
      SUBROUTINE COORD_NUM_SIMPLE(GBLEN,Gbin,dxd,dxd_br,XX)

      IMPLICIT NONE
        INTEGER(KIND=int16),INTENT(IN) :: GBLEN
        INTEGER(KIND=int16),INTENT(INOUT) :: Gbin(GBLEN,2)
        REAL(KIND=r8),INTENT(IN) :: dxd,dxd_br,XX(natms2)

!local variables
       INTEGER(KIND=int16) :: i,ix,iy,jx,jy,coord
       REAL(KIND=r8) :: dxsd,dxsd2,dxd_brn,dxd_brn2,xix,yix,yjx,xij,yij,rij2
       REAL(KIND=r8) :: dxbd,dxbd2

       frzymin=0.0_r8
       dxsd=dxd; dxsd2=dxsd*dxsd
       dxd_brn=dxd_br + 0.2_r8; dxd_brn2=dxd_brn*dxd_brn
       dxbd=dxd+dxss+0.2; dxbd2=dxbd*dxbd

       DO i=1,GBLEN

          coord=0_int16; brnatms=0_int16; br1natms=0_int16
          ix=Gbin(i,1); iy=ix+natms; ngbharr(1)=ix
          xix=XX(ix); yix=XX(iy)
          DO jx=1,natms
             jy=jx+natms
             yjx=XX(jy)
             xij=xix-XX(jx); xij=xij-boxl*ANINT(xij*rboxl)
             yij=yix-yjx
             rij2=xij*xij+yij*yij
             IF(jx .NE. ix) THEN
               IF (rij2 .LE. dxsd2) THEN
                   coord=coord+1_int16
                   ngbharr(coord+1)=jx
               ELSE IF(rij2 .LE. dxd_brn2 ) THEN
                   IF(rij2 <= dxbd2) THEN
                      br1natms=br1natms+1_int16
                      ngbh1arr_br(br1natms)=jx
                   END IF
                   brnatms=brnatms+1
                   ngbharr_br(brnatms)=jx
                   IF(yjx <= frzymin)frzymin=yjx
               END IF
             END IF
           END DO
           Gbin(i,2)=coord
          
      END DO
      intnatms=coord+1

      END SUBROUTINE COORD_NUM_SIMPLE
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
      SUBROUTINE FULL_POT_ENE(XX,Uene)

      IMPLICIT NONE
        REAL(KIND=r8),INTENT(IN)  :: XX(natms2)
        REAL(KIND=r8),INTENT(OUT) :: Uene
!local variables
       INTEGER(KIND=int16) :: ix,iy,jx,jy
       REAL(KIND=r8) :: xi,yi,si,ei,xij,yij,eij,rrij2,rrij6,sij,srij
       REAL(KIND=r8) :: sij6,dvdr,rijc,scrij
                        
       Uene=0.0_r8
       DO ix=1,natms-1
          iy=ix+natms; xi=XX(ix); yi=XX(iy) 
          si=sig(ix); ei=eps(ix)
          DO jx=ix+1,natms
             jy=jx+natms

             xij=xi-XX(jx)
             xij=xij-boxl*ANINT(xij*rboxl)
             yij=yi-XX(jy)
            rijc=rcut-SQRT(xij*xij+yij*yij)
           IF(rijc >= 0.0_r8) THEN

             rrij2=1.0_r8/(xij*xij + yij*yij)
             rrij6=rrij2*rrij2*rrij2
 
             sij=0.5_r8*(si+sig(jx))
             srij=sij**6 *rrij6

             eij=4.0_r8*srij*sqrt(ei*eps(jx))
             Uene=Uene + eij*(srij-1.0_r8)

      !       sij6=sij**6
      !       srij=sij6 *rrij6
      !       scrij=sij*rrcut6

      !       eij=4.0_r8*sqrt(ei*eps(jx))
      !       dvdr=-6.0_r8*eij*rrcut2*scrij*(2.0_r8*scrij-1.0_r8)
      !       Uene=Uene + eij*(srij*(srij-1.0_r8)-scrij*(scrij-1.0_r8))+dvdr*rijc
            END IF
          END DO
       END DO
 
       END SUBROUTINE FULL_POT_ENE
!!!!!!!!!!!!!!!!!!
      SUBROUTINE POT_DIST(ix,XX,Udist,Uene)

      IMPLICIT NONE
        INTEGER(KIND=int16) :: ix
        REAL(KIND=r8),INTENT(IN)  :: XX(natms2)
        REAL(KIND=r8),INTENT(OUT) :: Uene,Udist
!local variables
       INTEGER(KIND=int16) :: j,iy,jx,jy
       REAL(KIND=r8) :: xi,yi,si,ei,xij,yij,eij,rrij2,rrij6,sij,srij,seij
       Uene=0.0_r8
       Udist=0.0_r8

          iy=ix+natms; xi=XX(ix); yi=XX(iy) 
          si=sig(ix); ei=eps(ix) 

       DO j=2,intnatms
          jx=ngbharr(j)
          jy=jx+natms
          xij=xi-XX(jx)
          xij=xij-boxl*ANINT(xij*rboxl)
          yij=yi-XX(jy)
          rrij2=1.0_r8/(xij*xij + yij*yij)
          rrij6=rrij2*rrij2*rrij2

          sij=0.5_r8*(si+sig(jx))
          srij=sij**6 *rrij6
          seij=sqrt(ei*eps(jx))
          eij=4.0_r8*srij*seij*(srij-1.0_r8)
          Uene=Uene+eij
          Udist=Udist+eij+seij
       END DO

       DO j=1,brnatms
          jx=ngbharr_br(j)
          jy=jx+natms

          xij=xi-XX(jx)
          xij=xij-boxl*ANINT(xij*rboxl)
          yij=yi-XX(jy)
          rrij2=1.0_r8/(xij*xij + yij*yij)
          rrij6=rrij2*rrij2*rrij2
          sij=0.5_r8*(si+sig(jx))
          srij=sij**6 *rrij6

          eij=4.0_r8*srij*sqrt(ei*eps(jx))
          Uene=Uene + eij*(srij-1.0_r8)
        END DO

      END SUBROUTINE POT_DIST
!!!!!!!!!!!!!!!!!!
      SUBROUTINE NEW_COORD(ix,iix,ndiv,Mx,My,dxs,XX)
      IMPLICIT NONE
        INTEGER(KIND=int16),INTENT(IN) :: ix,iix,ndiv,Mx,My
        REAL(KIND=r8),INTENT(IN)  :: dxs
        REAL(KIND=r8),INTENT(INOUT)  :: XX(natms2)
!local variables
       INTEGER(KIND=int16) :: i,j,iy,iiy,jx,jy,mm,Abin(1,2),TestSpots(ndiv)
       REAL(KIND=r8) :: xi,yi,si,ei,xij,yij,eij,rrij2,rrij6,sij,srij
       REAL(KIND=r8) :: xc,yc,xix,yix,Uene,dxd,dxd_br,xitemp,yitemp


       EneArr=0.0_r8
  
       HEADT=HEAD; LISTT=LIST 
      
       iy=ix+natms; iiy=iix+natms
       xc=XX(iix); yc=XX(iiy)
       si=sig(ix); ei=eps(ix); xix=xc; yix=yc
       xitemp=xx(ix); yitemp=xx(iy)
       DO i=1,ndiv

          Uene=0.0_r8;
          xi=xc+NewR(i,1); yi=yc+NewR(i,2)
          xi=xi-boxl*ANINT(xi*rboxl)
          XX(ix)=xi; XX(iy)=yi
          CALL LINKS_ONE(ix,Mx,My,xx,XX)
!          CALL LINKS(Mx,My,XX)
          Abin(1,1)=ix; Abin(1,2)=0_int16
          dxd=REAL(5_int16,KIND=r8)*dxs; dxd_br=dxd+dxs
          CALL COORD_NUM(1_int16,Abin,Mx,My,dxd,dxd_br,xx)
          HEAD=HEADT; LIST=LISTT
           
          DO j=2,intnatms
             jx=ngbharr(j)
             jy=jx+natms

             xij=xi-XX(jx);
             xij=xij-boxl*ANINT(xij*rboxl)
             yij=yi-XX(jy)
             rrij2=1.0_r8/(xij*xij + yij*yij)
             rrij6=rrij2*rrij2*rrij2

             sij=0.5_r8*(si+sig(jx))
             srij=sij**6 *rrij6

             eij=4.0_r8*srij*sqrt(ei*eps(jx))
             Uene=Uene + eij*(srij-1.0_r8)
          END DO
          DO j=1,brnatms
            jx=ngbharr_br(j)
            jy=jx+natms

            xij=xi-XX(jx)
            xij=xij-boxl*ANINT(xij*rboxl)
            yij=yi-XX(jy)
            rrij2=1.0_r8/(xij*xij + yij*yij)
            rrij6=rrij2*rrij2*rrij2

            sij=0.5_r8*(si+sig(jx))
            srij=sij**6 *rrij6

            eij=4.0_r8*srij*sqrt(ei*eps(jx))
            Uene=Uene + eij*(srij-1.0_r8)
          END DO
          EneArr(i)=Uene
       END DO
       XX(ix)=xitemp; XX(iy)=yitemp
      END SUBROUTINE NEW_COORD
!!!!!!!!!!!!!!!!!!
      SUBROUTINE POT_1_NGBH(XX,Uene)

      IMPLICIT NONE
        REAL(KIND=r8),INTENT(IN)  :: XX(natms2)
        REAL(KIND=r8),INTENT(OUT) :: Uene
!local variables
       INTEGER(KIND=int16) :: i,j,ix,iy,jx,jy
       REAL(KIND=r8) :: xi,yi,si,ei,xij,yij,eij,rrij2,rrij6,sij,srij

       Uene=0.0_r8

       DO i=1,intnatms
          ix=ngbharr(i); iy=ix+natms; xi=XX(ix); yi=XX(iy) 
          si=sig(ix); ei=eps(ix)
          DO j=i+1,intnatms
             jx=ngbharr(j)
             jy=jx+natms

             xij=xi-XX(jx);
             xij=xij-boxl*ANINT(xij*rboxl)
             yij=yi-XX(jy)
             rrij2=1.0_r8/(xij*xij + yij*yij)
             rrij6=rrij2*rrij2*rrij2

             sij=0.5_r8*(si+sig(jx))
             srij=sij**6 *rrij6

             eij=4.0_r8*srij*sqrt(ei*eps(jx))
             Uene=Uene + eij*(srij-1.0_r8)
          END DO

          DO j=1,brnatms
            jx=ngbharr_br(j)
            jy=jx+natms

            xij=xi-XX(jx)
            xij=xij-boxl*ANINT(xij*rboxl)
            yij=yi-XX(jy)
            rrij2=1.0_r8/(xij*xij + yij*yij)
            rrij6=rrij2*rrij2*rrij2

            sij=0.5_r8*(si+sig(jx))
            srij=sij**6 *rrij6

            eij=4.0_r8*srij*sqrt(ei*eps(jx))
            Uene=Uene + eij*(srij-1.0_r8)
          END DO
       END DO


      END SUBROUTINE POT_1_NGBH
!!!!!!!!!!!!!!!!!!
      SUBROUTINE SPLINE(x,x1,x2,f1,f2,p1,p2,ff)
!  Program file name: spline.f90                                          !
!                                                                         !
!  <A9> Tao Pang 2006                                                        !
!                                                                         !
!  Last modified: June 23, 2006                                           !
!                                                                         !
!  (1) This F90 program is created for the book, "An Introduction to      !
!      Computational Physics, 2nd Edition," written by Tao Pang and       !
!      published by Cambridge University Press on January 19, 2006.       !
!                                                                         !
!  (2) No warranties, express or implied, are made for this program.      !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
     IMPLICIT NONE
     REAL(KIND=r8),INTENT(IN) :: x,x1,x2,f1,f2,p1,p2
     REAL(KIND=r8),INTENT(INOUT) :: ff

!Local variables
     REAL(KIND=r8) ::  DX, ALPHA, BETA, GAMMA, ETA

! Find the value of function f(x)
       DX = x2 - x1
       ALPHA = p2/(6.0_r8*DX)
       BETA = -P1/(6.0_r8*DX)
       GAMMA = f2/DX - DX*P2/6.0_r8
       ETA = DX*P1/6.0_r8 - f1/DX
       FF = ALPHA*(x-x1)*(x-x1)*(x-x1) &
          +BETA*(x-x2)*(x-x2)*(x-x2) &
          +GAMMA*(x-x1)+ETA*(x-x2)

      END SUBROUTINE SPLINE
!!!!!!!!!!!!!!!!!!
      SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
!  (1) This F90 program is created for the book, "An Introduction to      !
!      Computational Physics, 2nd Edition," written by Tao Pang and       !
!      published by Cambridge University Press on January 19, 2006.       !
!                                                                         !
!  N : number of data points
!  XI: x data points
!  FI: y data points
!  P2: 

  IMPLICIT NONE 
  INTEGER(KIND=int16) :: I
  INTEGER(KIND=int16), INTENT (IN) :: N
  REAL(KIND=r8), INTENT (IN), DIMENSION (N):: XI, FI
  REAL(KIND=r8), INTENT (OUT), DIMENSION (N):: P2
  REAL(KIND=r8), DIMENSION (N):: G, H
  REAL(KIND=r8), DIMENSION (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  DO I = 1, N-1
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO
!
! Evaluate the coefficient matrix elements
  DO I = 1, N-2
    D(I) = 2.0_r8*(H(I+1)+H(I))
    B(I) = 6.0_r8*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO
!
! Obtain the second-order derivatives
!
  CALL TRIDIAGONAL_LINEAR_EQ (N-2, D, C, C, B, G)
  P2(1) = 0.0_r8
  P2(N) = 0.0_r8
  DO I = 2, N-1 
    P2(I) = G(I-1)
  END DO
      END SUBROUTINE CUBIC_SPLINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER(KIND=int16), INTENT (IN) :: L
  INTEGER(KIND=int16) :: I
  REAL(KIND=r8), INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL(KIND=r8), INTENT (OUT), DIMENSION (L):: Z
  REAL(KIND=r8), DIMENSION (L):: Y, W
  REAL(KIND=r8), DIMENSION (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
      END SUBROUTINE TRIDIAGONAL_LINEAR_EQ
!!!!!!!!!!!!!!!!!!
      SUBROUTINE BINARY_SEARCH(value,a,found)
!Adapted from Programmer's guide to Fortran 90
! By Walter S. Brainerd, Charles H. Goldberg, Jeanne C. Adams

       implicit none
        real(kind=r8), intent(in) ::value,a(:)
        INTEGER(KIND=int16), intent(out) :: found

        INTEGER(KIND=int16) :: first,half,last, only

        first=1
        last=size(a)

        DO
          IF(first == last)exit
          half=(first+last)/2
          if(value <= a(half))then
            last=half
          else
            first = half+1
          end if
        END DO
        only=first
        if(value <= a(only)) then
          found=only
        else
          found=0
        end if
     END SUBROUTINE BINARY_SEARCH
!!!!!!!!!!!!!!!!!!
function binarySearch_I (a, value)
!
    INTEGER(KIND=int16)                  :: binarySearch_I
    real(kind=r8), intent(in), target :: a(:)
    real(kind=r8), intent(in)         :: value
    real(kind=r8), pointer            :: p(:)
    INTEGER(KIND=int16)                  :: mid, offset

    p => a
    binarySearch_I = 0
    offset = 0
    do while (size(p) > 0)
        mid = size(p)/2 + 1
        if ((p(mid) > value) .AND. (size(p) > 1) ) then
            p => p(:mid-1)
        else if ((p(mid) < value) .AND. (size(p) >1)) then
            offset = offset + mid
            p => p(mid+1:)
        else 
            binarySearch_I = offset + mid    ! SUCCESS!!
            return
        end if
    end do
end function binarySearch_I

!!!!!!!!!!!!!!!!!!
      SUBROUTINE TTIME(timebeg,timeend,totaltime)

      IMPLICIT NONE
!
! TTIME computes the time difference in seconds between
! the timestamps TIMEBEG and TIMEEND returned by the
! f90 intrinsic DATE_AND_TIME
!
      INTEGER(KIND=int16),DIMENSION(8),INTENT(INOUT) :: timebeg,timeend
      REAL(KIND=r8),INTENT(OUT) :: totaltime

! TIMEEND is modifed by borrowing in case each of its fields
! are not .GE. to the corresponding field in TIMEBEG (up to
! and including days)

      IF (timeend(8) .LT. timebeg(8)) THEN
          timeend(8)=timeend(8)+1000
          timeend(7)=timeend(7)-1
      END IF
      IF (timeend(7) .LT. timebeg(7)) THEN
          timeend(7)=timeend(7)+60
          timeend(6)=timeend(6)-1
      END IF
      IF (timeend(6) .LT. timebeg(6)) THEN
          timeend(6)=timeend(6)+60
          timeend(5)=timeend(5)-1
      END IF
      IF (timeend(5) .LT. timebeg(5)) THEN
          timeend(5)=timeend(5)+24
          timeend(3)=timeend(3)-1
      END IF

      totaltime=  REAL(timeend(8)-timebeg(8),KIND=r8) + &
            1000.0_r8*( REAL(timeend(7)-timebeg(7),KIND=r8) +&
              60.0_r8*( REAL(timeend(6)-timebeg(6),KIND=r8) +&
              60.0_r8*( REAL(timeend(5)-timebeg(5),KIND=r8) +&
              24.0_r8*( REAL(timeend(3)-timebeg(3),KIND=r8)))))
      totaltime=totaltime/1000.0_r8

      RETURN
      END SUBROUTINE TTIME

      function duni()

      implicit none

      logical new
      INTEGER(KIND=int16) ir,jr,i,j,k,l,m,ii,jj
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

      if(new)then

!     initial values of i,j,k must be in range 1 to 178 (not all 1)
!     initial value of l must be in range 0 to 168.

        i=12
        j=34
        k=56
        l=78
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

!     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
      return

  end function duni 

      subroutine gauss(natms,vxx,vyy)
!=============================================================
      implicit none

      INTEGER(KIND=int16) natms,i
      real(8) vxx,vyy,a1,a3,a5,a7,a9,randum,rrr,rr2
      
      dimension vxx(natms),vyy(natms)
      
      data a1,a3,a5/3.949846138d0,0.252408784d0,0.076542912d0/
      data a7,a9/0.008355968d0,0.029899776d0/

!     initialise random number generator

      randum=duni()

      do i=1,natms
        
        rrr=(duni()+duni()+duni()+duni()+duni()+duni()&
         +duni()+duni()+duni()+duni()+duni()+duni()&
         -6.d0)/4.d0
        rr2=rrr*rrr
        vxx(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
        
        rrr=(duni()+duni()+duni()+duni()+duni()+duni()&
         +duni()+duni()+duni()+duni()+duni()+duni() &
         -6.d0)/4.d0
        rr2=rrr*rrr
        vyy(i)=rrr*(a1+rr2*(a3+rr2*(a5+rr2*(a7+rr2*a9))))
        
      enddo
      return
      end subroutine gauss 


end module main 
