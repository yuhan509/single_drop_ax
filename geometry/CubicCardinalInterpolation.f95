SUBROUTINE CubicCardinalInterpolation_arrayo()
!!// ith element, jth node, kth abscissa point
!!// s1o is first filled in here
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER:: j,i,k
DOUBLE PRECISION :: s,ps1i,tmp,xr,xm

Do i = 1, N
     ps1i = ps1(i)
     xr = (ps1(i+1)-ps1i)*.5d0
     xm = (ps1(i+1)+ps1i)*.5d0    
     Do  k = 1, qm
          s = xm + xr * qx(k)
          s1o(k,i) = s
          tmp = s - ps1i
          
          Do j = 1, N+1    
               
               car_na1o(j,k,i) = acar_na(j,1,i)*tmp**3 + acar_na(j,2,i)*tmp**2 &
                + acar_na(j,3,i)*tmp
               IF (i == j) car_na1o(j,k,i) = car_na1o(j,k,i) + 1d0
               dcards_na1o(j,k,i) = 3d0*acar_na(j,1,i)*tmp**2 + 2d0*acar_na(j,2,i)*tmp &
                + acar_na(j,3,i)
               d2cards2_na1o(j,k,i) = 6d0*acar_na(j,1,i)*tmp + 2d0*acar_na(j,2,i)
               
               car_cl1o(j,k,i) = acar_cl(j,1,i)*tmp**3 + acar_cl(j,2,i)*tmp**2 &
                + acar_cl(j,3,i)*tmp
               IF (i == j) car_cl1o(j,k,i) = car_cl1o(j,k,i) + 1d0
               dcards_cl1o(j,k,i) = 3d0*acar_cl(j,1,i)*tmp**2 + 2d0*acar_cl(j,2,i)*tmp &
                + acar_cl(j,3,i)
               d2cards2_cl1o(j,k,i) = 6d0*acar_cl(j,1,i)*tmp + 2d0*acar_cl(j,2,i)

          ENDDO
     ENDDO
ENDDO
               
END SUBROUTINE
               
               
               
SUBROUTINE CubicCardinalInterpolation_special_arrayo()
!!// ith element, jth node, kth abscissa point
!!// s1o_s is first filled in here
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER:: j,i,k
DOUBLE PRECISION :: s,tmp,xr

Do i = 1, N

     xr = (ps1(i+1)-ps1(i))     
     Do  k = 1, 2*sqm
     
          IF (k .le. sqm) THEN
              s = xr * sqx(k) + ps1(i)
          ELSE
              s = ps1(i+1) - xr * sqx(k-sqm)
          ENDIF         
          s1o_s(k,i) = s
          tmp = s - ps1(i)
          
          Do j = 1, N+1    
               
               car_na1o_s(j,k,i) = acar_na(j,1,i)*tmp**3 + acar_na(j,2,i)*tmp**2 &
                + acar_na(j,3,i)*tmp
               IF (i == j) car_na1o_s(j,k,i) = car_na1o_s(j,k,i) + 1d0
               dcards_na1o_s(j,k,i) = 3d0*acar_na(j,1,i)*tmp**2 + 2d0*acar_na(j,2,i)*tmp &
                + acar_na(j,3,i)
               d2cards2_na1o_s(j,k,i) = 6d0*acar_na(j,1,i)*tmp + 2d0*acar_na(j,2,i)
               
               car_cl1o_s(j,k,i) = acar_cl(j,1,i)*tmp**3 + acar_cl(j,2,i)*tmp**2 &
                + acar_cl(j,3,i)*tmp
               IF (i == j) car_cl1o_s(j,k,i) = car_cl1o_s(j,k,i) + 1d0
               dcards_cl1o_s(j,k,i) = 3d0*acar_cl(j,1,i)*tmp**2 + 2d0*acar_cl(j,2,i)*tmp &
                + acar_cl(j,3,i)
               d2cards2_cl1o_s(j,k,i) = 6d0*acar_cl(j,1,i)*tmp + 2d0*acar_cl(j,2,i)

          ENDDO
     ENDDO
ENDDO
               
END SUBROUTINE
