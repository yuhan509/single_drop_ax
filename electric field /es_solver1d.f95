SUBROUTINE es_Solver1d()
! use N,r,z,ar,az,acar,dphidn,PI
use mod_elasbody
IMPLICIT NONE
INTEGER,PARAMETER :: NRHS=1,N1=N+1,N2= N+2!,N3=2*N+3,N4=2*N+4
INTEGER :: INFO
INTEGER :: i0,j,ipiv(N2) ! 2*N+2 unknowns
DOUBLE PRECISION :: AM(N2,N2),BM(N2,NRHS)
DOUBLE PRECISION :: phi1!,phi2!,start,finish
DOUBLE PRECISION :: c11,c12,c3!,c21,c22
 !c1_ for i0,j on same drop; c2_ for i0,j on diff drops

AM = 0d0
Do i0=1,N1
    !AM(i0,N+2) = 2d0*PI ! careful with the sign    
    Do j=1,N1
        CALL es_dcoef_gax(i0,j,c11)
        !write(*,*) 'after AM i0 = 1  j = 1 gax'
        AM(i0,j) = c11
!        AM(i0+N1,j+N1) = c11 !mirror  
!        CALL es_dcoef2d_gax(i0,j,c21)
!        AM(i0,j+N1) = c21
!        AM(i0+N1,j) = c21 !mirror
    ENDDO
    
    CALL es_dcoef_dgdnax(i0,c12)
    AM(i0,N2) = 2d0*PI - c12  !!coef of phi1
!    AM(i0+N1,N4) = 2d0*PI - c12  !!mirror coef of phi2
!    CALL es_dcoef2d_dgdnax(i0,c22)
!    AM(i0,N4) = -c22          !!coef of phi2
!    AM(i0+N1,N3) = -c22       !!mirror coef of phi1

    BM(i0,NRHS) = -4d0*PI*z1(i0)!!phi_inf
!    BM(i0+N1,NRHS) = 4d0*PI*z2(i0)
END DO

!At the interface net charge = 0,then
BM(N2,NRHS) = 0d0 !! 0 net charge
!BM(N4,NRHS) = 0d0 !! 0 net charge   
Do j=1,N1
   call es_ich_coef1(j,c3) 
   AM(N2,j) = c3
!   AM(N4,j+N1) = c3
ENDDO

!Do i=1,N+2
!   Do j=1,N+2
!   write(60,*) AM(i,j)!"AM,i,j",AM(i,j),i,j
!   ENDDO
!ENDDO

!DO i=1,N+2
!write(61,*) BM(i,NRHS)!"BM",BM
!enddo
!write(*,*) 'after AM'
CALL dgesv(N2,NRHS,AM,N2,IPIV,BM,N2,INFO)

!IF(INFO==0) THEN
DO j=1,N1
   dphidn1(j) = BM(j,NRHS)
!   dphidn2(j) = BM(j+N1,NRHS)
ENDDO
!write(62,*) 
!phi = BM(N+2,NRHS)
phi1 = BM(N2,NRHS) 
!phi2 = BM(N4,NRHS)
write(15,*) phi1!,',',phi2
print*, 'phi1= ',phi1,'INFO = ',INFO!,'phi2= ',phi2,'INFO = ',INFO

!ELSE
!WRITE (*,*) "Error using es_solver!!, INFO=",INFO
!END IFes_dcoef.f95

END SUBROUTINE

