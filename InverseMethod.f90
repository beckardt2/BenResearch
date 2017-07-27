Program Power_Method

IMPLICIT NONE
Interface 
	Subroutine DominantEigen(N,dim,x,lamda)
		integer, intent(in)::dim
		Double Precision, intent(in), dimension(dim,dim) :: N
		Double Precision, intent(out), dimension(dim)::x
		Real, intent(out)::lamda
	End Subroutine DominantEigen

	Subroutine LUfactor(N,dim,L,U)
		integer, intent(in)::dim
                Double Precision, intent(in), dimension(dim,dim) :: N
		Double Precision, intent(out), dimension(dim,dim) :: L, U
	End Subroutine LUfactor
	Subroutine LeastDominant(L,U,dim,vec,val)
		integer, intent(in)::dim
		Double Precision, intent(in), dimension(dim,dim) :: L, U
		Double Precision, intent(out), dimension(dim) :: vec
		Real, intent(out) :: val
	End Subroutine LeastDominant
End Interface

double precision, allocatable :: A(:,:), lower(:,:), upper(:,:), inverse(:,:)
integer :: i,j,ndim, err
Double Precision, allocatable,  DIMENSION(:) :: vector
real:: value
Integer:: INFO
double precision, allocatable:: VR(:,:),WR(:),WI(:),WORK(:),VL(:),pos(:)


open(15, file="HESSIAN.txt")

read(15,*) ndim
allocate(A(ndim,ndim), lower(ndim,ndim), upper(ndim,ndim), inverse(ndim,ndim),VL(ndim),WR(ndim),WI(ndim),VR(ndim,ndim)&
,WORK(4*ndim),pos(ndim),  Stat=err)
IF (err/=0) Then
        STOP
End if
read(15,*) ((A(i,j), j=1,ndim), i=1,ndim)
close (15)

allocate(vector(ndim), Stat=err)
IF (err/=0) Then
        STOP
End if

!Call DominantEigen(A,ndim,vector,value)

!print*, "The dominant Eigenvector of the given matrix is", vector, " with corresponding dominant Eigenvalue", value

!Call LUfactor(A,ndim,lower,upper)
!print*,lower
!print*,upper

!open(unit=1, file='L_U2.txt', ACTION="write", STATUS="replace")
!        do i=1,ndim
!            write(1,*) real( lower(i,:) )
!        end do
!        write(1,*) 'upper'
!        do i=1,ndim
!            write(1,*) real( upper(i,:) )
!        end do
!        close(1)



!Call LeastDominant (lower,upper,ndim,vector,value)

!print*, "weakest Eigvalue", value
!print*, vector
print*,"LAPACK says:"
Call DGEEV('N','V',ndim,A,ndim,WR,WI,VL,ndim,VR,ndim,WORK,4*ndim,INFO)
!Print*, VR
print*, WR

print*," Dominant eigenvector"
print*,"Dominant eigenvalue", MAXVAL(WR)
open(unit=2, file='DOM.txt', ACTION="write", STATUS="replace")
write(2,*) real( VR(:,MAXLOC(WR)))
close(2)

Do i=1,ndim
	pos(i)=abs(WR(i))
End Do
print*," weakest eigenvalue ", WR(MINLOC(pos))

open(unit=3, file='WEAK.txt', ACTION="write", STATUS="replace")
write(3,*) real( VR(:,MINLOC(WR)))
close(3)


End Program Power_Method

SUBROUTINE DominantEigen(N,dim,x,lamda)

Implicit none
integer, intent(in)::dim
Double Precision, intent(in), dimension(dim,dim) :: N
Double Precision, intent(out), dimension(dim)::x
Real, intent(out)::lamda
Real, dimension(dim)::xi,xf,xc
real::i,c

xi(:)=1
xf=matmul(N,xi)

x=xf/maxval(xf)
xc=xi
!the purpose of xc is to compare to x to see how much it changed
Do While (abs(dot_product(x-xc,x-xc))>0.0000001)
	xc=x
	xi=xf
	xf=matmul(N,xi)
        x=xf/xf(dim)
End Do

lamda= dot_product(x,matmul(N,x))/dot_product(x,x)

End SUBROUTINE DominantEigen

Subroutine LUfactor(N,dim,L,U)

Implicit none
integer, intent(in)::dim
Double Precision, intent(in), dimension(dim,dim) :: N
Double Precision, intent(out), dimension(dim,dim) :: L,U
integer:: i,j,k

Do i=1,dim-1,1
	L(i,i+1:)=0
	U(i+1:,i)=0
End Do

Do i=1,dim,1
	L(i,1)=N(i,1)
	U(1,i)=N(1,i)/N(1,1)
	U(i,i)=1
End Do

Do j=2,dim-1,1
	Do i=j,dim,1
		L(i,j)= N(i,j)
		Do k=1,j-1,1
			L(i,j)=L(i,j)-L(i,k)*U(k,j)
		End Do
	End Do
	Do k=j+1,dim,1
		U(j,k)=N(j,k)
		Do i=1,j-1,1
			U(j,k)=U(j,k)-L(j,i)*U(i,k)
		end Do
		U(j,k)=U(j,k)/L(j,j)
	end Do
end do
Do i=1,dim,1
	L(i,i)=N(i,i)
	Do j=1,i-1
		L(i,i)=L(i,i)-L(i,j)*U(j,i)
	End Do
End Do

End subroutine LUfactor

subroutine LeastDominant(L,U,dim,vec,val)

implicit none
integer, intent(in)::dim
Double Precision, intent(in), dimension(dim,dim) :: L, U
Double Precision, intent(out), dimension(dim) :: vec
real, intent(out) :: val
Double Precision, dimension(dim)::xi,xf,y,d,posx,posy
integer:: i,j,k
real:: mew, v,large

xi(:)=0
xf(:)=1
Do While (abs(dot_product(xf-xi,xf-xi))>0.0000001)
	xi=xf
	Do j=1,dim,1
		d(j)=xi(j)
		Do k=1,j-1,1
			d(j)=d(j)-L(j,k)*d(k)
		End Do
		d(j)=d(j)/L(j,j)
	End Do
	Do j=0,dim-1,1
		y(dim-j)=d(dim-j)
		Do k=1,j,1
			y(dim-j)=y(dim-j)-y(dim-k+1)*U(dim-j,dim-k+1)
		End Do
	End Do
	
	IF (abs(maxval(y))>abs(minval(y))) THEN
                mew=maxval(y)
        ELSE
                 mew=minval(y)
        END IF

	
	
	v=1/mew
	xf=v*y
	Do i=1, dim, 1
                posx(i)=abs(xf(i))
        End Do
        large=maxval(posx)
	Do i=1, dim, 1
                xf(i)=xf(i)/large
        End Do
End Do
vec=xf
val=v

End Subroutine LeastDominant 
		




