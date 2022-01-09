SUBROUTINE DilaRateMeasure1d(istep)
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER,INTENT(IN) :: istep
INTEGER :: i,k,j
DOUBLE PRECISION :: s,rs,trs,tzs,gu22,gu11,car_na,dcards_na,dcards_cl
DOUBLE PRECISION :: cj,ck,DilaRate,t22_dv,t11_dv
DOUBLE PRECISION :: cj_1,x11,t11_sv,t22_sv

IF (mod(istep,10) == 1) THEN

     write(19,*) 'istep = ', istep,' s, dilatation rate, t11_dv part, t22_dv part'
     write(29,*) 'istep = ', istep,' s, shear rate related, t11_sv part, t22_sv part'
     Do i = 1,N
          Do k = 1,qm
               s = s1o(k,i)
               rs = r1o(k,i)
               trs = tr1o(k,i)
               tzs = tz1o(k,i)
               gu22 = gu22o1(k,i)
               gu11 = gu11o1(k,i)
               DilaRate = 0d0
               x11 = 0d0
               
               Do j = 1, N+1    
                    
                    car_na = car_na1o(j,k,i)
                    dcards_na = dcards_na1o(j,k,i)
                    dcards_cl = dcards_cl1o(j,k,i)
                                   
                    cj = car_na/rs + gu22 * trs * dcards_na
                    ck = gu22 * tzs * dcards_cl
                    !! DilaRate = g^11 * S_11 + g^22 * S_22                 
                    DilaRate = DilaRate + cj * ur1(j) + ck * uz1(j)

                    cj_1 = car_na/rs - gu22 * trs * dcards_na
                    !! x11 = g^11 * S_11 - g^22 * S_22
                    x11 = x11 + cj_1 * ur1(j) - ck * uz1(j)
                  
               Enddo
               t11_dv = DilaRate * dilaViscosity * gu11
               t22_dv = DilaRate * dilaViscosity * gu22

               t11_sv = x11 * gu11
               t22_sv = -x11 * gu22
               write(19,300) s,',', DilaRate,',',t11_dv,',',t22_dv
               write(29,300) s,',', x11,',',t11_sv,',',t22_sv
          Enddo
     Enddo

ENDIF

100 FORMAT(ES25.15,A1,ES25.15)
101 FORMAT(I8,A1,ES25.15,A1,ES25.15)
200 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)
300 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15,A1,ES25.15)
END SUBROUTINE
