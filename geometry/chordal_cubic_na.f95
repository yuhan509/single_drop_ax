SUBROUTINE chordal_cubics_na(N,ps,x,ax)
implicit NONE
INTEGER,INTENT(IN) :: N
DOUBLE PRECISION,INTENT(IN),DIMENSION(N+1) :: x,ps
DOUBLE PRECISION,DIMENSION(N+1) :: bx1
DOUBLE PRECISION,INTENT(OUT),DIMENSION(3,N) :: ax
DOUBLE PRECISION,DIMENSION(N+1) :: h,a,b,c,yx
!a,b,c contains diagonal entries of a (n)*(n) matrix
!a is the diagonal entry 
!b and c are the upper and lower diagonal respectively
INTEGER :: i

Do i=1,N
h(i)=ps(i+1)-ps(i)
ENDDO

Do i=1,N-1
a(i+1)=(h(i)+h(i+1))*2d0/3d0
b(i+1)=h(i+1)/3d0
c(i)=h(i)/3d0
yx(i+1)=(x(i+2)-x(i+1))/h(i+1)-(x(i+1)-x(i))/h(i)
ENDDO

!!
a(1)=1d0
yx(1)=0d0
a(N+1)=1d0
yx(N+1)=0d0
!! natural cubic spline

Call Thomas_algorithm(N+1,a,b,c,yx,bx1)

Do i=1,N
ax(1,i)=(bx1(i+1)-bx1(i))/(3d0*h(i))
ax(2,i)=bx1(i)   !!discard the last entry of bx1
ax(3,i)=(x(i+1)-x(i))/h(i)-h(i)/3d0*(bx1(i+1)+2d0*bx1(i))
ENDDO

END SUBROUTINE
