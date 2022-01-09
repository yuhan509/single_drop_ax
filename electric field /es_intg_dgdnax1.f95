DOUBLE PRECISION FUNCTION es_intg_dgdnax1o(i0,i,k,s)
! use r1,z1
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0,i,k!,j
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION :: Fm,Em,a,b,m1,zd,dgdnax!,car
DOUBLE PRECISION :: rs,nrs,zs,nzs,dlds

!    rs = es_r1o(k,i) 
!    zs = es_z1o(k,i) 
!    dlds = es_dlds1o(k,i)
!    nrs = es_nr1o(k,i)
!    nzs = es_nz1o(k,i)
     rs= r1o(k,i) 
     zs= z1o(k,i) 
     dlds= dlds1o(k,i) 
     nrs = nr1o(k,i)
     nzs = nz1o(k,i)
     
zd = zs-z1(i0)
a = r1(i0)*r1(i0)+rs*rs+zd*zd
b = 2d0*rs*r1(i0)
m1 = (a-b)/(a+b)
CALL Ellip_int_dnons(i0,i,s,a,b,m1,Fm,Em)
!write(*,*) 'dgdnax1','i0',i0,'i',i,'s',s,'Fm',Fm,'Em',Em
!call getcar(s,j,car)

dgdnax=2d0*nrs/(rs*dsqrt(a+b))*(Em-Fm)&
 -4d0*(nrs*(rs-r1(i0))+nzs*zd)/((a-b)*dsqrt(a+b))*Em
es_intg_dgdnax1o=rs*dlds*dgdnax!*car

END FUNCTION


DOUBLE PRECISION FUNCTION es_intg_sdgdnax1(i0,i,s)
! use r1,z1
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0,i!,j
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION :: sFm,a,b,zd,sdgdnax!,car
DOUBLE PRECISION :: rs,nrs,zs,nzs,dlds

Call geom_cubs_21(i,s,rs,nrs,zs,nzs,dlds)
Call Ellip_int_dsing(i0,sFm)
!Call getcar(s,j,car)

zd = zs-z1(i0)
a = r1(i0)*r1(i0)+rs*rs+zd*zd
b = 2d0*rs*r1(i0)
sdgdnax = 2d0*nrs/(rs*dsqrt(a+b))*(-sFm)
es_intg_sdgdnax1=rs*dlds*sdgdnax!*car

END FUNCTION
