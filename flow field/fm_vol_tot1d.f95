subroutine fm_vol_tot1d(istep)
!use N
!also evaluate the instantaneous centroid pos along z axis
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: istep
INTEGER :: i
DOUBLE PRECISION :: ai,bi,ss,totvol,totarea!,centroidz
DOUBLE PRECISION,EXTERNAL :: vol_elem1d,area_elem1d!,cen_elem1d

totvol = 0d0
totarea = 0d0
!centroidz = 0d0
Do i = 1,N
   ai = ps1(i)
   bi = ps1(i+1)
   call vol_gq_cubs(vol_elem1d,i,ai,bi,ss)
   totvol = totvol + ss
   call vol_gq_cubs(area_elem1d,i,ai,bi,ss)
   totarea = totarea + ss
!   call vol_gq_cubs(cen_elem1d,i,ai,bi,ss)
!   centroidz = centroidz + ss
ENDDO
write(53,*) istep,',',totvol*PI,',',totarea*2d0*PI
!write(54,*) centroidz/totvol
END SUBROUTINE



SUBROUTINE vol_gq_cubs(func,i,ai,bi,ss)
IMPLICIT NONE
!! a and b are the lower and upper borders, a=i,b=i+1
DOUBLE PRECISION,EXTERNAL :: func
INTEGER,INTENT(IN) :: i
DOUBLE PRECISION,INTENT(IN) :: ai,bi !integral range
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
w = (/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,.0666713443d0/)
x = (/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,.9739065285d0/)
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
    ss = ss + w(k)*(func(i,xm+dx)+func(i,xm-dx))
ENDDO
ss = ss*xr
!write(*,*) "regular_gauss xr",xr
END SUBROUTINE


DOUBLE PRECISION FUNCTION vol_elem1d(i,s)
!use N,r,z,ar,az
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
INTEGER,INTENT(IN) :: i
double PRECISION :: temp,rs,dzds

temp = s-ps1(i)
rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)

vol_elem1d = rs*rs*dzds
END FUNCTION


DOUBLE PRECISION FUNCTION area_elem1d(i,s)
!use N,r,z,ar,az
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
INTEGER,INTENT(IN) :: i
double PRECISION :: temp,rs,dzds,drds

temp = s-ps1(i)
rs = ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
dzds = 3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)
drds = 3d0*ar1(1,i)*temp**2+2d0*ar1(2,i)*temp+ar1(3,i)

area_elem1d = rs*dsqrt(dzds**2+drds**2)

END FUNCTION


DOUBLE PRECISION FUNCTION cen_elem1d(i,s)
!use N,r,z,ar,az
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i
DOUBLE PRECISION,INTENT(IN) :: s
double PRECISION :: temp,rs,zs,dzds

temp = s-ps1(i)
rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
zs=az1(1,i)*temp**3+az1(2,i)*temp**2+az1(3,i)*temp+z1(i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)

cen_elem1d = rs*rs*zs*dzds
END FUNCTION
