SUBROUTINE initialContour

USE mod_elasbody
IMPLICIT NONE

CALL sphere_nodes_rz()
CALL arc_para()
!!!! cubic cardinal function and derivatives
!!// first gives all abscissa 
CALL CubicCardinalInterpolation_arrayo()
CALL CubicCardinalInterpolation_special_arrayo()

!!!! rest reference metrics
call geo1_rest_metrics_arrayo_cubs_gq()
call geom_special_rest_metrics_arrayo_cubs_gq()
!!!! 


100 FORMAT(ES25.15,A1,ES25.15)
101 FORMAT(I8,A1,ES25.15,A1,ES25.15)
200 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)

END SUBROUTINE



SUBROUTINE recorder(istep)

USE mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: istep
INTEGER :: j

    write(40,*) 'istep=',istep
    write(22,*) 'istep=',istep,'uz, ur'
    write(23,*) 'istep=',istep
!    write(72,*) 'istep=',istep
!    write(16,*) 'istep=',istep,'regular gauss abscissa'
!    write(17,*) 'istep=',istep,'special gauss abscissa'
    
    do j = 1,N+1
        write(23,100) dphidn1(j)
        write(21,100) z1(j),',',r1(j)
        write(22,100) uz1(j),',',ur1(j)  
        write(40,200) ps1(j),',',z1(j),',',r1(j)         
    enddo
    write(*,*) istep,',',(z1(N+1)-r1(N/2+1))/(z1(N+1)+r1(N/2+1))
    write(13,*) istep,',',(z1(N+1)-r1(N/2+1))/(z1(N+1)+r1(N/2+1))    
    
    CALL fm_vol_tot1d(istep)
    
    CALL DilaRateMeasure1d(istep)
    
100 FORMAT(ES25.15,A1,ES25.15)
101 FORMAT(I8,A1,ES25.15,A1,ES25.15)
200 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)

END SUBROUTINE
