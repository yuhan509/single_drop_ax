SUBROUTINE Ellip_int_dsing(i0,sFm) 
! use N,SN
use mod_elasbody
!use mod_singbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0
DOUBLE PRECISION,INTENT(OUT) :: sFm
DOUBLE PRECISION :: b0


IF(i0==1 .or. i0==N+1) THEN
    sFm = 0d0
ELSE
    b0 = .5d0
    sFm = -2d0*b0 ! sFm: factor before ln(s-s0)
ENDIF
         
END SUBROUTINE


!SUBROUTINE Ellip_int_dsing_sb(i0,sFm) 
!! use N,SN
!use mod_singbody
!IMPLICIT NONE
!INTEGER,INTENT(IN) :: i0
!DOUBLE PRECISION,INTENT(OUT) :: sFm
!DOUBLE PRECISION :: b0
!
!IF(i0==1 .or. i0==SN+1) THEN
!    sFm = 0d0
!ELSE
!    b0 = .5d0
!    sFm = -2d0*b0 ! sFm: factor before ln(s-s0)
!ENDIF
!         
!END SUBROUTINE
