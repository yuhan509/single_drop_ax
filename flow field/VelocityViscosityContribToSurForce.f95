
!!// called only after contour shifted
SUBROUTINE coefForUnodej_arrayo()
!!// ith element, jth node, kth abscissa point
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER :: j,i,k
DOUBLE PRECISION :: cj,ck,cl,cm,cpr,cqr,cpz,cqz
DOUBLE PRECISION :: rs,trs,tzs,gu22,dgu22ds,d2rds2,d2zds2
DOUBLE PRECISION :: car_na,dcards_na,d2cards2_na,car_cl,dcards_cl,d2cards2_cl
DOUBLE PRECISION :: nrs,nzs,partofcovarintderiv,meancurv


Do i = 1, N
     Do  k = 1, qm
          rs = r1o(k,i)
          trs = tr1o(k,i)
          tzs = tz1o(k,i)
          d2rds2 = d2rds2o1(k,i)
          d2zds2 = d2zds2o1(k,i)
          gu22 = gu22o1(k,i)
          dgu22ds = dgu22dso1(k,i)
          
          nrs = nr1o(k,i)
          nzs = nz1o(k,i)
          partofcovarintderiv = partofcovarintderiv1o(k,i)
          meancurv = meancurv1o(k,i)
          
          
          Do j = 1, N+1    
               
               car_na = car_na1o(j,k,i)
               dcards_na = dcards_na1o(j,k,i)
               d2cards2_na = d2cards2_na1o(j,k,i)
               
               car_cl = car_cl1o(j,k,i)
               dcards_cl = dcards_cl1o(j,k,i)
               d2cards2_cl = d2cards2_cl1o(j,k,i)
               
               cj = car_na/rs + gu22 * trs * dcards_na
               ck = gu22 * tzs * dcards_cl
               cl = -trs/(rs**2)*car_na + (d2rds2*gu22 + trs*dgu22ds)*dcards_na &
                    + trs*gu22*d2cards2_na
               cm = (d2zds2*gu22 + tzs*dgu22ds)*dcards_cl + tzs*gu22*d2cards2_cl
               
!! sign of cp,cq
               cpr = -dilaViscosity * (trs * partofcovarintderiv + nrs * meancurv)
               cqr = -dilaViscosity * trs * gu22
               
               cpz = -dilaViscosity * (tzs * partofcovarintderiv + nzs * meancurv)
               cqz = -dilaViscosity * tzs * gu22
               
               
               delForce_r_coefUrj1o(j,k,i) = cpr * cj + cqr * cl
               delForce_r_coefUzj1o(j,k,i) = cpr * ck + cqr * cm
               delForce_z_coefUrj1o(j,k,i) = cpz * cj + cqz * cl
               delForce_z_coefUzj1o(j,k,i) = cpz * ck + cqz * cm
               
!               write(85,*) 'j=',j,'i=',i,'k=',k,'dcards_na=',dcards_na,'dcards_cl',dcards_cl
!               write(84,*) 'j=',j,'i=',i,'k=',k,'car_na=',car_na,'car_cl',car_cl
          ENDDO
     ENDDO
ENDDO
               
END SUBROUTINE
               
               
               
SUBROUTINE coefForUnodej_special_arrayo()
!!// ith element, jth node, kth abscissa point
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER :: j,i,k
DOUBLE PRECISION :: cj,ck,cl,cm,cpr,cqr,cpz,cqz
DOUBLE PRECISION :: rs,trs,tzs,gu22,dgu22ds,d2rds2,d2zds2
DOUBLE PRECISION :: car_na,dcards_na,d2cards2_na,car_cl,dcards_cl,d2cards2_cl
DOUBLE PRECISION :: nrs,nzs,partofcovarintderiv,meancurv

Do i = 1, N
     Do  k = 1, 2*sqm           
          rs = r1o_s(k,i)
          trs = tr1o_s(k,i)
          tzs = tz1o_s(k,i)
          d2rds2 = d2rds2o1_s(k,i)
          d2zds2 = d2zds2o1_s(k,i)
          gu22 = gu22o1_s(k,i)
          dgu22ds = dgu22dso1_s(k,i)
          
          nrs = nr1o_s(k,i)
          nzs = nz1o_s(k,i)
          partofcovarintderiv = partofcovarintderiv1o_s(k,i)
          meancurv = meancurv1o_s(k,i)

          Do j = 1, N+1
               
               car_na = car_na1o_s(j,k,i)
               dcards_na = dcards_na1o_s(j,k,i)
               d2cards2_na = d2cards2_na1o_s(j,k,i)
               
               car_cl = car_cl1o_s(j,k,i)
               dcards_cl = dcards_cl1o_s(j,k,i)
               d2cards2_cl = d2cards2_cl1o_s(j,k,i)
               
               cj = car_na/rs + gu22 * trs * dcards_na
               ck = gu22 * tzs * dcards_cl
               cl = -trs/(rs**2)*car_na + (d2rds2*gu22 + trs*dgu22ds)*dcards_na &
                    + trs*gu22*d2cards2_na
               cm = (d2zds2*gu22 + tzs*dgu22ds)*dcards_cl + tzs*gu22*d2cards2_cl
               
!! sign of cp,cq
               cpr = -dilaViscosity * (trs * partofcovarintderiv + nrs * meancurv)
               cqr = -dilaViscosity * trs * gu22
               
               cpz = -dilaViscosity * (tzs * partofcovarintderiv + nzs * meancurv)
               cqz = -dilaViscosity * tzs * gu22
                              
               delForce_r_coefUrj1o_s(j,k,i) = cpr * cj + cqr * cl
               delForce_r_coefUzj1o_s(j,k,i) = cpr * ck + cqr * cm
               delForce_z_coefUrj1o_s(j,k,i) = cpz * cj + cqz * cl
               delForce_z_coefUzj1o_s(j,k,i) = cpz * ck + cqz * cm
              

          ENDDO
     ENDDO
ENDDO

END SUBROUTINE
