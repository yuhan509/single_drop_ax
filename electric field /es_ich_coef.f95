SUBROUTINE es_ich_coef1(j,c3) 
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j
DOUBLE PRECISION,INTENT(OUT) :: c3
DOUBLE PRECISION :: ai,bi,ss !,temp
DOUBLE PRECISION,EXTERNAL :: intg_cubs01
INTEGER :: i

c3 = 0d0
Do i=1,N ! ith elemental integral 
    ai = ps1(i)
    bi = ps1(i+1)
    CALL card_gq_cubs(intg_cubs01,i,ai,bi,j,ss)
    c3 = c3 + ss
ENDDO
END SUBROUTINE


DOUBLE PRECISION FUNCTION intg_cubs01(i,s,j)
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,j !i for ith element; j for cardinal contribution on jth node
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION :: car
DOUBLE PRECISION :: rs,dlds

Call geom_cubs_01(i,s,rs,dlds)
Call getcar_cl(i,s,j,car)

intg_cubs01=car*rs*dlds
END FUNCTION


SUBROUTINE card_gq_cubs(func,i,ai,bi,j,ss)
!! a and b are the lower and upper borders, a=i,b=i+1
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
DOUBLE PRECISION,INTENT(IN) :: ai,bi !integral range
INTEGER,INTENT(IN) :: i,j
DOUBLE PRECISION,INTENT(OUT) :: ss
INTEGER :: k
INTEGER,PARAMETER :: m=5
DOUBLE PRECISION :: dx,xr,xm,w(m),x(m) 
!x = (/.095012509837637440185d0,.281603550779258913230d0,&
!.458016777657227586342d0,.617876244402643748447d0,&
!.755404408355003033895d0,.865631202387831743880d0,&
!.944575023073232576078d0,.989400934991649932596d0/)
!w = (/.189450610455068496285d0,.182603415044923588867d0,&
!.169156519395002538189d0,.149595988816576732081d0,&
!.124628971255533872052d0,.095158511682492784810d0,&
!.062253523938647892863d0,.027152459411754094852d0/)
!
w = (/.2955242247d0,.2692667193d0,.2190863625d0,&
.1494513491d0,.0666713443d0/)
x = (/.1488743389d0,.4333953941d0,.6794095682d0,&
.8650633666d0,.9739065285d0/)
!
!x=(/.339981d0, .861136d0 /)
!w=(/.652145d0, .347855d0 /)
!
!x=(/.577350269189626d0/)
!w=(/1d0/)

ss = 0d0
xr = (bi-ai)*.5d0
xm = (bi+ai)*.5d0
Do k=1,m
    dx = xr * x(k)
    ss = ss + w(k)*(func(i,xm+dx,j)+func(i,xm-dx,j))
ENDDO
ss = ss*xr
!write(*,*) "regular_gauss xr",xr
END SUBROUTINE

