SUBROUTINE Thomas_algorithm(n,a,b,c,y,x)
!a,b,c contains diagonal entries of a (n)*(n) matrix
!a is the diagonal entry
!b and c are the upper and lower diagonal respectively
INTEGER,INTENT(IN) :: n
DOUBLE PRECISION,INTENT(IN),DIMENSION(n) :: a,y
DOUBLE PRECISION,INTENT(IN),DIMENSION(n-1) :: b,c
DOUBLE PRECISION,DIMENSION(n) :: s
DOUBLE PRECISION,DIMENSION(n-1) :: d
DOUBLE PRECISION,INTENT(OUT),DIMENSION(n) :: x
INTEGER :: i

d(1)=b(1)/a(1)
s(1)=y(1)/a(1)
Do i=1,n-1
d(i+1)=1/(a(i+1)-c(i)*d(i))*b(i+1)
s(i+1)=1/(a(i+1)-c(i)*d(i))*(y(i+1)-c(i)*s(i))
ENDDO

x(n)=s(n)
Do i=1,n-1
x(n-i)=s(n-i)-d(n-i)*x(n-i+1)
ENDDO

END SUBROUTINE






