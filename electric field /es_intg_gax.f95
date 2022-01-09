DOUBLE PRECISION FUNCTION es_intg_gax1o(j,i0,i,k,s)
! use r1,z1
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i0,i,k
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION :: car,Fm,Em,m1,a,b,zd,gax
DOUBLE PRECISION :: rs,zs,dlds

!    rs = es_r1o(k,i) 
!    zs = es_z1o(k,i) 
!    dlds = es_dlds1o(k,i)
    rs = r1o(k,i) 
    zs = z1o(k,i) 
    dlds = dlds1o(k,i) 

zd = zs-z1(i0)
a = r1(i0)*r1(i0)+rs*rs+zd*zd
b = 2d0*rs*r1(i0)
m1 = (a-b)/(a+b)
CALL Ellip_int_dnons(i0,i,s,a,b,m1,Fm,Em)
call getcar_cl(i,s,j,car)

gax = 4d0*Fm/dsqrt(a+b)
es_intg_gax1o=car*rs*dlds*gax

END FUNCTION


DOUBLE PRECISION FUNCTION es_intg_sgax1(j,i0,i,s)
! use r1,z1
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i0,i
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION :: car,sFm,a,b,zd,sgax
DOUBLE PRECISION :: rs,zs,dlds

Call geom_cubs_11(i,s,rs,zs,dlds)
Call Ellip_int_dsing(i0,sFm)
Call getcar_cl(i,s,j,car)

zd = zs-z1(i0)
a = r1(i0)*r1(i0)+rs*rs+zd*zd
b = 2d0*rs*r1(i0)
sgax = 4d0*sFm/dsqrt(a+b)
es_intg_sgax1=car*rs*dlds*sgax

END FUNCTION
