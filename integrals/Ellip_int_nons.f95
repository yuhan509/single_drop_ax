SUBROUTINE Ellip_int_dnons(i0,i,s,a,b,m1,Fm,Em)
! calculate ellipitc integrals using polynomial approximation
! use N,ps1
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0,i
DOUBLE PRECISION,INTENT(IN) :: s,a,b,m1
DOUBLE PRECISION,INTENT(OUT) :: Fm,Em

DOUBLE PRECISION :: a0=1.38629436112d0,&
 a1=.09666344259d0, a2=.03590092383d0,&
 a3=.03742563713d0, a4=.01451196212d0,&
 b0=.5d0, b1=.12498593597d0, b2=.06880248576d0,&
 b3=.03328355346d0, b4=.00441787012d0,&
 c1=.44325141463d0, c2=.06260601220d0,& 
 c3=.04757383546d0, c4=.01736506451d0,&
 d1=.24998368310d0, d2=.09200180037d0,&
 d3=.04069697526d0, d4=.00526449639d0

DOUBLE PRECISION :: temp

! i0 represents kernel node

Em = (1d0+c1*m1+c2*m1**2d0+c3*m1**3d0+c4*m1**4d0)& 
-(d1*m1+d2*m1**2d0+d3*m1**3d0+d4*m1**4d0)*DLOG(m1)

!IF(i0==1 .or. i0==N+1 .or. i0==SN+1) THEN
IF(i0==1 .or. i0==N+1) THEN
    Fm = a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0
ELSE IF(i0.eq.i .or. i0.eq.i+1) THEN
    temp = (s-ps1(i))/(ps1(i+1)-ps1(i))+i-i0  !!somehow fix it
    Fm = (a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0)&
    -(b1*m1+b2*m1**2d0+b3*m1**3d0+b4*m1**4d0)*DLOG(m1)&
    +b0*DLOG(a+b)-b0*DLOG((a-b)/temp**2)
ELSE
    Fm = (a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0)& 
    -(b0+b1*m1+b2*m1**2d0+b3*m1**3d0+b4*m1**4d0)*DLOG(m1)    
ENDIF

END SUBROUTINE



!SUBROUTINE Ellip_int_dnons2d(m1,Fm,Em)
!! calculate ellipitc integrals using polynomial approximation
!! use N
!use mod_elasbody
!IMPLICIT NONE
!DOUBLE PRECISION,INTENT(IN) :: m1
!DOUBLE PRECISION,INTENT(OUT) :: Fm,Em
!
!DOUBLE PRECISION :: a0=1.38629436112d0,&
! a1=.09666344259d0, a2=.03590092383d0,&
! a3=.03742563713d0, a4=.01451196212d0,&
! b0=.5d0, b1=.12498593597d0, b2=.06880248576d0,&
! b3=.03328355346d0, b4=.00441787012d0,&
! c1=.44325141463d0, c2=.06260601220d0,& 
! c3=.04757383546d0, c4=.01736506451d0,&
! d1=.24998368310d0, d2=.09200180037d0,&
! d3=.04069697526d0, d4=.00526449639d0
!
!Em = (1d0+c1*m1+c2*m1**2d0+c3*m1**3d0+c4*m1**4d0)& 
!-(d1*m1+d2*m1**2d0+d3*m1**3d0+d4*m1**4d0)*DLOG(m1)
!
!Fm = (a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0)& 
!-(b0+b1*m1+b2*m1**2d0+b3*m1**3d0+b4*m1**4d0)*DLOG(m1)    
!
!END SUBROUTINE


!SUBROUTINE Ellip_int_dnons_sb(i0,i,s,a,b,m1,Fm,Em)
!! calculate ellipitc integrals using polynomial approximation
!! use N,SN,psb
!use mod_singbody
!IMPLICIT NONE
!INTEGER,INTENT(IN) :: i0,i
!DOUBLE PRECISION,INTENT(IN) :: s,a,b,m1
!DOUBLE PRECISION,INTENT(OUT) :: Fm,Em
!
!DOUBLE PRECISION :: a0=1.38629436112d0,&
! a1=.09666344259d0, a2=.03590092383d0,&
! a3=.03742563713d0, a4=.01451196212d0,&
! b0=.5d0, b1=.12498593597d0, b2=.06880248576d0,&
! b3=.03328355346d0, b4=.00441787012d0,&
! c1=.44325141463d0, c2=.06260601220d0,& 
! c3=.04757383546d0, c4=.01736506451d0,&
! d1=.24998368310d0, d2=.09200180037d0,&
! d3=.04069697526d0, d4=.00526449639d0
!
!DOUBLE PRECISION :: temp
!
!! i0 represents kernel node
!
!Em = (1d0+c1*m1+c2*m1**2d0+c3*m1**3d0+c4*m1**4d0)& 
!-(d1*m1+d2*m1**2d0+d3*m1**3d0+d4*m1**4d0)*DLOG(m1)
!
!IF(i0==1 .or. i0==SN+1) THEN
!    Fm = a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0
!ELSE IF(i0.eq.i .or. i0.eq.i+1) THEN
!    temp = (s-psb(i))/(psb(i+1)-psb(i))+i-i0  !!somehow fix it
!    Fm = (a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0)&
!    -(b1*m1+b2*m1**2d0+b3*m1**3d0+b4*m1**4d0)*DLOG(m1)&
!    +b0*DLOG(a+b)-b0*DLOG((a-b)/temp**2)
!ELSE
!    Fm = (a0+a1*m1+a2*m1**2d0+a3*m1**3d0+a4*m1**4d0)& 
!    -(b0+b1*m1+b2*m1**2d0+b3*m1**3d0+b4*m1**4d0)*DLOG(m1)    
!ENDIF
!
!END SUBROUTINE
