SUBROUTINE solveVelocity(istep)

!USE mod_gaussqxw 
USE mod_elasbody

IMPLICIT NONE
INTEGER,INTENT(IN) :: istep

CALL geo1_metrics_arrayo_cubs_gq(istep)
CALL geo1_special_metrics_arrayo_cubs_sgq()

CALL coefForUnodej_arrayo()
CALL coefForUnodej_special_arrayo()
!write(*,*) 'after coefForUnodej'
CALL es_Solver1d()
!write(*,*) 'after es_Solver'
CALL for1_arrayo_cubs_gq()
CALL for1_special_arrayo_cubs_gq()
!write(*,*) 'after for1_arrayo'
CALL USolver()


100 FORMAT(ES25.15,A1,ES25.15)
101 FORMAT(I8,A1,ES25.15,A1,ES25.15)
200 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)

END SUBROUTINE
