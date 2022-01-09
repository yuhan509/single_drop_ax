SUBROUTINE USolver()
! use N,ur1,uz1,ur2,uz2,PI
! viscosity ratio = 1.0
use mod_elasbody
IMPLICIT NONE
INTEGER,PARAMETER :: NRHS = 1,UN1= N+1,UN2 = 2*N+2
INTEGER :: INFO
INTEGER :: i0,j,ipiv(UN2) ! 2*N+2 unknowns
DOUBLE PRECISION :: UAM(UN2,UN2),UBM(UN2,NRHS)
DOUBLE PRECISION :: curr1,curz1,cuzr1,cuzz1,fr1,fz1

UAM = 0d0
Do i0=1,N+1
    !AM(i0,N+2) = 2d0*PI ! careful with the sign    
    Call UGaxisymmetry_arrayo_cubs_gq(i0)
    Call UGaxisymmetry_special_arrayo_cubs_sgq(i0)
    
    Do j=1,N+1
        CALL Ucoef_gax(i0,j,curr1,curz1,cuzr1,cuzz1)
        IF(j==i0) Then
            curr1 = curr1 + 8d0*pi
            cuzz1 = cuzz1 + 8d0*pi
        ENDIF
        UAM(i0,j)=curr1
        UAM(i0,j+N+1)=curz1
        UAM(i0+N+1,j)=cuzr1
        UAM(i0+N+1,j+N+1)=cuzz1
    ENDDO
    !!Terms didn't involve Uz, Ur
    CALL Fcoef_gax(i0,fr1,fz1)  
    UBM(i0,NRHS) = - fr1 
    UBM(i0+N+1,NRHS) = - fz1 
END DO

!Do i=1,N+2
!   Do j=1,N+2
!   write(60,*) AM(i,j)!"AM,i,j",AM(i,j),i,j
!   ENDDO
!ENDDO

!DO i=1,N+2
!write(61,*) BM(i,NRHS)!"BM",BM
!enddo
CALL dgesv(UN2,NRHS,UAM,UN2,IPIV,UBM,UN2,INFO)

!IF(INFO==0) THEN
DO j=1,N+1
   ur1(j) = UBM(j,NRHS)
   uz1(j) = UBM(j+N+1,NRHS)
   
ENDDO


END SUBROUTINE
