SUBROUTINE es_sgq_gax(func,i,ai,bi,j,i0,ss)
!!! Integrand form: f(x)*log(x)=sum of f(x)
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
DOUBLE PRECISION,INTENT(IN) :: ai,bi !integral range
INTEGER,INTENT(IN) :: i0,j,i
DOUBLE PRECISION,INTENT(OUT) :: ss
INTEGER :: k
INTEGER,PARAMETER :: m=4
DOUBLE PRECISION :: w(m),x(m),s,dx 

w = (/-0.383464d0,-0.386875d0,-0.190435d0,-0.039225d0/)
x = (/.041448d0,.245275d0,.556165d0,.848982d0/)

!w = (/-0.718539d0, -0.281461d0/)
!x = (/0.112009d0, 0.602277d0/)

!write(*,*) 'INT(ai)',i
ss=0d0
Do k=1,m
    dx = (bi-ai)*x(k)
    IF(i==i0) THEN
    s = dx+ai
    ELSE IF(i+1==i0) THEN
    s = -dx+bi
    END IF 
    !write(*,*) 'sgq s=',s
    ss = ss + w(k)*func(j,i0,i,s)
    !write(*,*) '*sgq ss=',ss
    !write(*,*) 'integrand2s',integrand_cubs2s(j,i0,s,N,r,z,ar,az,acar)
END DO
ss = ss*(bi-ai)

END SUBROUTINE



SUBROUTINE es_sgq_dgdnax(func,i,ai,bi,i0,ss)
!!! Integrand form: f(x)*log(x)=sum of f(x)
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
DOUBLE PRECISION,INTENT(IN) :: ai,bi !integral range
INTEGER,INTENT(IN) :: i0,i
DOUBLE PRECISION,INTENT(OUT) :: ss
INTEGER :: k
INTEGER,PARAMETER :: m=4
DOUBLE PRECISION :: w(m),x(m),s,dx 

w = (/-0.383464d0,-0.386875d0,-0.190435d0,-0.039225d0/)
x = (/.041448d0,.245275d0,.556165d0,.848982d0/)

!w = (/-0.718539d0, -0.281461d0/)
!x = (/0.112009d0, 0.602277d0/)

!write(*,*) 'INT(ai)',i
ss=0d0

Do k=1,m
    dx = (bi-ai)*x(k)
    IF(i==i0) THEN
    s = dx+ai
    ELSE IF(i+1==i0) THEN
    s = -dx+bi
    END IF 
    !write(*,*) 'sgq s=',s
    ss = ss + w(k)*func(i0,i,s)
    !write(*,*) '*sgq ss=',ss
    !write(*,*) 'integrand2s',integrand_cubs2s(j,i0,s,N,r,z,ar,az,acar)
END DO
ss = ss*(bi-ai)

END SUBROUTINE

