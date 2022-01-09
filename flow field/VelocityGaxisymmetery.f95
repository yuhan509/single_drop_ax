!//To be called in solver when i0 changing
SUBROUTINE UGaxisymmetry_arrayo_cubs_gq(i0)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0

INTEGER :: i,k
DOUBLE PRECISION :: zd,zd2,km,a,b,c,d,f,g,h,m,m1!,car,temp
DOUBLE PRECISION :: s,rs,zs
DOUBLE PRECISION :: Fm,Em,Mrz,Mrr,Mzz,Mzr
DOUBLE PRECISION :: I30,I10!,I32

Do i = 1,N
    IF (i0 == N+1 .and. i == 1) THEN
         write(*,*) '*'
    ENDIF
     Do k = 1,qm
          s = s1o(k,i)
          rs = r1o(k,i) 
          zs = z1o(k,i) 
          
          zd = zs - z1(i0)
          zd2 = zd**2d0
          d = r1(i0)**2 - rs**2
          f = r1(i0)**2 + rs**2
          
          a = f + zd2
          b = 2d0*rs*r1(i0)
          c = dsqrt(rs/r1(i0))
          g = a - b
          h = dsqrt(a+b)
          
          m = 2d0*b/(a+b) 
          m1 = 1d0-m
          km = dsqrt(m)
          
          Call Ellip_int_dnons(i0,i,s,a,b,m1,Fm,Em)
          !print*, 'Fm',Fm,'Em',Em
          ! Integrate Stokeslet over azimuthal angle 0 to 2*PI
          ! Integrate component of G
          IF(i0==1 .or. i0==N+1) THEN
              I10 = 2d0/h*PI
              I30 = 2d0/h**3d0*PI
              !I32 = 1d0/h**3d0*PI
              Mzz = rs*(I10+zd2*I30)
              Mzr = rs**2d0*zd*I30
              Mrz = 0d0
              Mrr = 0d0
          Else
              Mzz = 2d0*km*c*(Fm+zd2*Em/g)
              Mzr = km/(c*r1(i0))*zd*(Fm-(d+zd2)*Em/g)
              Mrz = -km*zd/r1(i0)*c*(Fm+(d-zd2)*Em/g)
              Mrr = 2d0*km/b*c*((a+zd2)*Fm-(2d0*zd2**2d0+3d0*zd2*f+&
              d**2d0)*Em/g)
          Endif
          
          Mzz1o(k,i)=Mzz
          Mzr1o(k,i)=Mzr
          Mrz1o(k,i)=Mrz
          Mrr1o(k,i)=Mrr
     ENDDO
ENDDO
!print*, 1,'Mrz',Mrz,'Mrr',Mrr,'Mzz',Mzz,'Mzr',Mzr

END SUBROUTINE



SUBROUTINE UGaxisymmetry_special_arrayo_cubs_sgq(i0)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0
INTEGER :: i,ks,firstks
DOUBLE PRECISION :: zd,zd2,k,a,b,c,d,f,g,h,m!,car,temp
DOUBLE PRECISION :: rs,zs
DOUBLE PRECISION :: sFm,sMzr,sMrr,sMrz,sMzz
DOUBLE PRECISION :: I30,I10!,I32


Do i = 1,N
     IF (i == i0 .or. i + 1 == i0) THEN
          IF (i == i0) THEN
              firstks = 1
          ELSE
              firstks = 1 + sqm
          ENDIF
          Do ks = firstks,firstks + sqm - 1
               
               rs = r1o_s(ks,i) 
               zs = z1o_s(ks,i) 
               
               zd = zs - z1(i0)
               zd2 = zd**2d0
               d = r1(i0)**2 - rs**2
               f = r1(i0)**2 + rs**2
               
               a = f + zd2
               b = 2d0*rs*r1(i0)
               c = dsqrt(rs/r1(i0))
               g = a - b
               h = dsqrt(a+b)
               
               m = 2d0*b/(a+b) 
               k = dsqrt(m)
               
               Call Ellip_int_dsing(i0,sFm) 
               
               ! Integrate Stokeslet over azimuthal angle 0 to 2*PI
               ! Integrate component of G
               IF(i0==1 .or. i0==N+1) THEN
                   I10 = 2d0/h*PI
                   I30 = 2d0/h**3d0*PI
                   !I32 = 1d0/h**3d0*PI
                   sMzz = 0d0
                   sMzr = 0d0
                   sMrz = 0d0
                   sMrr = 0d0
               Else
                   sMzz = 2d0*k*c*sFm
                   sMzr = k/(c*r1(i0))*zd*sFm
                   sMrz = -k*zd/r1(i0)*c*sFm
                   sMrr = 2d0*k/b*c*(a+zd2)*sFm
               ENDIF
               Mzz1o_s(ks,i) = sMzz
               Mzr1o_s(ks,i) = sMzr
               Mrz1o_s(ks,i) = sMrz
               Mrr1o_s(ks,i) = sMrr
          ENDDO
     ENDIF
ENDDO

END SUBROUTINE
