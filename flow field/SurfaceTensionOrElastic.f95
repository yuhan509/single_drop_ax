SUBROUTINE SurfaceTensionDependOnCoverage(js, djsds, gu11, gu22, dgu22ds,&
b11, b22, c112, c211, c222, normal_stress, tangent_stress) 
use mod_elasbody

DOUBLE PRECISION, INTENT(IN) :: js, djsds, gu11, gu22, dgu22ds, b11, b22,&
c112, c211, c222
DOUBLE PRECISION, INTENT(OUT) :: normal_stress, tangent_stress
DOUBLE PRECISION :: gamma, dt22ds, t11, t22

! Surface tension as function of surface coverage (in terms of js)
        gamma = 1d0 + surftenMAX * DLOG(1d0 - eta/js)
        dt22ds = gu22 * surftenMax / (1d0 - eta/js) * eta / (js * js) * djsds &
                + dgu22ds * gamma
!!        gamma = 1d0 + surftenMAX * DLOG(1d0 - eta)
!!        dt22ds = dgu22ds * gamma

!!// here only the contribution of surface tension is included
!!// because it is independent of velocity
        t11 = gamma*gu11
        t22 = gamma*gu22
        
        normal_stress = b11 * t11 + b22 * t22
        tangent_stress = dt22ds + (c112 + c222 + c222) * t22 + c211 * t11
        
END SUBROUTINE



SUBROUTINE SurfaceTensionConstant(gu11, gu22, dgu22ds,&
b11, b22, c112, c211, c222, normal_stress, tangent_stress) 
use mod_elasbody

DOUBLE PRECISION, INTENT(IN) :: gu11, gu22, dgu22ds, b11, b22,&
c112, c211, c222
DOUBLE PRECISION, INTENT(OUT) :: normal_stress, tangent_stress
DOUBLE PRECISION :: gamma, dt22ds, t11, t22

! Surface tension as function of surface coverage (in terms of js)
!!        gamma = 1d0 
        gamma = 1d0 + surftenMAX * DLOG(1d0 - eta)
        dt22ds = dgu22ds * gamma

!!// here only the contribution of surface tension is included
!!// because it is independent of velocity
        t11 = gamma*gu11
        t22 = gamma*gu22
        
        normal_stress = b11 * t11 + b22 * t22
        tangent_stress = dt22ds + (c112 + c222 + c222) * t22 + c211 * t11
        
END SUBROUTINE


SUBROUTINE ElasticStress(js, djsds, gu11, gu22, dgu22ds, ogu11_dum, ogu22_dum,& 
odgu22ds_dum, b11, b22, c112, c211, c222, normal_stress, tangent_stress) 

use mod_elasbody

DOUBLE PRECISION, INTENT(IN) :: js, djsds, gu11, gu22, dgu22ds,&
ogu11_dum, ogu22_dum, odgu22ds_dum, b11, b22, c112, c211, c222
DOUBLE PRECISION, INTENT(OUT) :: normal_stress, tangent_stress
DOUBLE PRECISION :: dt22ds, t11, t22

        t11 = gs * (1d0 / js * ogu11_dum - 1d0 / js**3 * gu11)
        t22 = gs * (1d0 / js * ogu22_dum - 1d0 / js**3 * gu22)
        dt22ds = gs * (-ogu22_dum / js**2 * djsds + odgu22ds_dum / js &
                + 3d0 / js**4 *djsds * gu22 - dgu22ds / js**3)

        normal_stress = b11 * t11 + b22 * t22
        tangent_stress = dt22ds + (c112 + c222 + c222) * t22 + c211 * t11
        
END SUBROUTINE
