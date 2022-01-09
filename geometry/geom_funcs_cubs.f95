SUBROUTINE geom_cubs_01(i,s,rs,dlds)
! use r1,z1,ar1,az1
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION,INTENT(OUT) :: rs,dlds
DOUBLE PRECISION :: temp,drds,dzds
INTEGER,INTENT(IN) :: i

temp = s-ps1(i)

rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
!first order derivatives
drds=3d0*ar1(1,i)*temp**2+2d0*ar1(2,i)*temp+ar1(3,i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)
dlds = dsqrt(drds*drds+dzds*dzds)
END SUBROUTINE 


SUBROUTINE geom_cubs_11(i,s,rs,zs,dlds)
! use r1,z1,ar1,az1
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION,INTENT(OUT) :: rs,zs,dlds
DOUBLE PRECISION :: temp,drds,dzds
INTEGER,INTENT(IN) :: i

temp = s-ps1(i)

rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
zs=az1(1,i)*temp**3+az1(2,i)*temp**2+az1(3,i)*temp+z1(i)
!first order derivatives
drds=3d0*ar1(1,i)*temp**2+2d0*ar1(2,i)*temp+ar1(3,i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)
dlds = dsqrt(drds**2d0+dzds**2d0)
END SUBROUTINE 



SUBROUTINE geom_cubs_21(i,s,rs,nrs,zs,nzs,dlds)
! use r1,z1,ar1,az1
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION,INTENT(OUT) :: rs,nrs,zs,nzs,dlds
DOUBLE PRECISION :: temp,drds,dzds
INTEGER,INTENT(IN) :: i

temp = s-ps1(i)

rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
zs=az1(1,i)*temp**3+az1(2,i)*temp**2+az1(3,i)*temp+z1(i)
!first order derivatives
drds=3d0*ar1(1,i)*temp**2+2d0*ar1(2,i)*temp+ar1(3,i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)
dlds = dsqrt(drds**2d0+dzds**2d0)

!normal component,pointing outward
nrs= dzds/dlds
nzs= -drds/dlds
!print*, 'nrs',nrs,'nzs',nzs
END SUBROUTINE 


SUBROUTINE geom_cubs_31(i,s,rs,nrs,zs,nzs,dlds,div_n)
! use r1,z1,ar1,az1
! only used in fm part
use mod_elasbody
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN) :: s
DOUBLE PRECISION,INTENT(OUT) :: rs,nrs,zs,nzs,dlds,div_n
DOUBLE PRECISION :: temp,temp1,&
dzds,drds,d2zds,d2rds,drdz,d2rdz
INTEGER,INTENT(IN) :: i

temp = s-ps1(i)

rs=ar1(1,i)*temp**3+ar1(2,i)*temp**2+ar1(3,i)*temp+r1(i)
zs=az1(1,i)*temp**3+az1(2,i)*temp**2+az1(3,i)*temp+z1(i)
!first order derivatives
drds=3d0*ar1(1,i)*temp**2+2d0*ar1(2,i)*temp+ar1(3,i)
dzds=3d0*az1(1,i)*temp**2+2d0*az1(2,i)*temp+az1(3,i)
dlds = dsqrt(drds**2d0+dzds**2d0)
!second order derivatives
d2rds = 6d0*ar1(1,i)*temp+2d0*ar1(2,i)
d2zds = 6d0*az1(1,i)*temp+2d0*az1(2,i)
!normal component,pointing outward
nrs= dzds/dlds
nzs= -drds/dlds
!print*, 'nrs',nrs,'nzs',nzs

!For axisymetric curvature
drdz = drds/dzds
 !second derivative of r in z
d2rdz = d2rds/(dzds*dzds)-d2zds*drds/(dzds**3)

temp1 = dsqrt(1d0+drdz*drdz)
! div of unit normal vector
! equal to twice the mean curvature
div_n = 1d0/(rs*temp1)-1d0/temp1**3*d2rdz

END SUBROUTINE 


!
!SUBROUTINE geom_cubs_12d(i,s,rs,zs,dlds)
!! use r2,z2,ar2,az2
!use mod_elasbody
!IMPLICIT NONE
!DOUBLE PRECISION,INTENT(IN) :: s
!DOUBLE PRECISION,INTENT(OUT) :: rs,zs,dlds
!DOUBLE PRECISION :: temp,drds,dzds
!INTEGER,INTENT(IN) :: i
!
!temp = s-ps1(i)
!
!rs=ar2(1,i)*temp**3+ar2(2,i)*temp**2+ar2(3,i)*temp+r2(i)
!zs=az2(1,i)*temp**3+az2(2,i)*temp**2+az2(3,i)*temp+z2(i)
!!first order derivatives
!drds=3d0*ar2(1,i)*temp**2+2d0*ar2(2,i)*temp+ar2(3,i)
!dzds=3d0*az2(1,i)*temp**2+2d0*az2(2,i)*temp+az2(3,i)
!dlds = dsqrt(drds**2d0+dzds**2d0)
!END SUBROUTINE 
!
!
!
!SUBROUTINE geom_cubs_22d(i,s,rs,nrs,zs,nzs,dlds)
!! use r2,z2,ar2,az2
!use mod_elasbody
!IMPLICIT NONE
!DOUBLE PRECISION,INTENT(IN) :: s
!DOUBLE PRECISION,INTENT(OUT) :: rs,nrs,zs,nzs,dlds
!DOUBLE PRECISION :: temp,drds,dzds
!INTEGER,INTENT(IN) :: i
!
!temp = s-ps1(i)
!
!rs=ar2(1,i)*temp**3+ar2(2,i)*temp**2+ar2(3,i)*temp+r2(i)
!zs=az2(1,i)*temp**3+az2(2,i)*temp**2+az2(3,i)*temp+z2(i)
!!first order derivatives
!drds=3d0*ar2(1,i)*temp**2+2d0*ar2(2,i)*temp+ar2(3,i)
!dzds=3d0*az2(1,i)*temp**2+2d0*az2(2,i)*temp+az2(3,i)
!dlds = dsqrt(drds**2d0+dzds**2d0)
!
!!normal component,pointing outward
!nrs= -dzds/dlds
!nzs= drds/dlds
!!print*, 'nrs',nrs,'nzs',nzs
!END SUBROUTINE 
!
!
!SUBROUTINE geom_cubs_32d(i,s,rs,nrs,zs,nzs,dlds,div_n)
!! use r2,z2,ar2,az2
!! only used in fm part
!use mod_elasbody
!IMPLICIT NONE
!DOUBLE PRECISION,INTENT(IN) :: s
!DOUBLE PRECISION,INTENT(OUT) :: rs,nrs,zs,nzs,dlds,div_n
!DOUBLE PRECISION :: temp,temp1,&
!dzds,drds,d2zds,d2rds,drdz,d2rdz
!INTEGER,INTENT(IN) :: i
!
!temp = s-ps1(i)
!
!rs=ar2(1,i)*temp**3+ar2(2,i)*temp**2+ar2(3,i)*temp+r2(i)
!zs=az2(1,i)*temp**3+az2(2,i)*temp**2+az2(3,i)*temp+z2(i)
!!first order derivatives
!drds=3d0*ar2(1,i)*temp**2+2d0*ar2(2,i)*temp+ar2(3,i)
!dzds=3d0*az2(1,i)*temp**2+2d0*az2(2,i)*temp+az2(3,i)
!dlds = dsqrt(drds**2d0+dzds**2d0)
!!second order derivatives
!d2rds = 6d0*ar2(1,i)*temp+2d0*ar2(2,i)
!d2zds = 6d0*az2(1,i)*temp+2d0*az2(2,i)
!!normal component,pointing outward
!nrs= -dzds/dlds
!nzs= drds/dlds
!!print*, 'nrs',nrs,'nzs',nzs
!
!!For axisymetric curvature
!drdz = drds/dzds
! !second derivative of r in z
!d2rdz = d2rds/(dzds*dzds)-d2zds*drds/(dzds**3)
!
!temp1 = dsqrt(1d0+drdz*drdz)
!! div of unit normal vector
!! equal to twice the mean curvature
!div_n = 1d0/(rs*temp1)-1d0/temp1**3*d2rdz
!
!END SUBROUTINE 
