SUBROUTINE geo1_special_metrics_arrayo_cubs_sgq()
! use dphidn1,N
use mod_elasbody 
use mod_gaussqxw
implicit none

integer :: i,k !!ith element !!kth special gq abscissa
DOUBLE PRECISION :: s,ps1i,rs,zs,dlds
DOUBLE PRECISION :: nrs,trs,nzs,tzs
DOUBLE PRECISION :: temp,dzds,drds,d2zds2,d2rds2
DOUBLE PRECISION :: ar11i,ar12i,ar13i,az11i,az12i,az13i,r1i,z1i
!! derivatives
!DOUBLE PRECISION :: drdl,dzdl,d2lds,d2rdl,d2zdl 
!! metric tensors for deformed state,g12 and g21 are zero
!! u means upper index, contravariant; d means lower index, covariant
DOUBLE PRECISION :: gu11, gu22, gd11, gd22
DOUBLE PRECISION :: dgu11ds, dgu22ds, dgd11ds, dgd22ds
!! christoffel symbol, c111, c122, c212 are zero
!! the first number is contravariant index and the last 2 numbers are covariant index
DOUBLE PRECISION :: c112, c211, c222
!! curvature tensor, b12 and b21 are zero
!! two numbers are covariant (lower) index
DOUBLE PRECISION :: b11, b22 
!! local extension ratio, js
DOUBLE PRECISION :: js, djsds
!! contravariant elastic stress resultants
!DOUBLE PRECISION :: t11, t22, dt22ds
!! surface tension, as a function of surface coverage
!DOUBLE PRECISION :: gamma
DOUBLE PRECISION :: normal_stress, tangent_stress


Do i = 1,N
    ar11i = ar1(1,i)
    ar12i = ar1(2,i) 
    ar13i = ar1(3,i)
    az11i = az1(1,i)
    az12i = az1(2,i)
    az13i = az1(3,i)
    r1i = r1(i)
    z1i = z1(i)
    ps1i = ps1(i) 
    Do k = 1,2*sqm
        s = s1o_s(k,i)
        temp = s - ps1i
                              
        rs=ar11i*temp**3+ar12i*temp**2+ar13i*temp+r1i
        zs=az11i*temp**3+az12i*temp**2+az13i*temp+z1i
        !first order derivatives
        drds=3d0*ar11i*temp**2+2d0*ar12i*temp+ar13i
        dzds=3d0*az11i*temp**2+2d0*az12i*temp+az13i
        dlds = dsqrt(drds**2+dzds**2)
!       !also unit tangent vector r,z component
!        drdl = drds / dlds 
!        dzdl = dzds / dlds
       
        !second order derivatives
        d2rds2 = 6d0*ar11i*temp+2d0*ar12i
        d2zds2 = 6d0*az11i*temp+2d0*az12i
             
!        d2lds = (drds * d2rds + dzds * d2zds) / dlds         
!        d2rdl = d2rds / (dlds*dlds) - drds * d2lds / (dlds**3)
!        d2zdl = d2zds / (dlds*dlds) - dzds * d2lds / (dlds**3)
 
        !normal vector r,z component,pointing outward
        nrs= dzds/dlds
        nzs= -drds/dlds
        ! a_2 vector r,z component        
        trs = drds
        tzs = dzds
                
        ! metric tensors
        gu11 = 1d0 / (rs*rs)
        gd11 = rs * rs
        gu22 = 1d0 / (drds**2 + dzds**2)        
        gd22 = drds**2 + dzds**2
        
        !dervivative of metric tensors
        dgd11ds = 2d0 * rs * drds
        dgd22ds = 2d0 * (drds * d2rds2 + dzds * d2zds2)
        dgu11ds = -2d0 * drds / rs**3
        dgu22ds = -2d0 * (drds * d2rds2 + dzds * d2zds2) / (drds**2 + dzds**2)**2 
                
        ! Christoffel symbols
        c112 = drds / rs
        c211 = -rs * drds / (drds**2 + dzds**2)
        c222 = (drds * d2rds2 + dzds * d2zds2) / (dzds**2 + drds**2)
        
        ! Curvature tensors
        b11 = - dzds * rs / dsqrt(dzds**2 + drds**2) 
        b22 = (dzds * d2rds2 - drds * d2zds2) / dsqrt(dzds**2 + drds**2)
        
        !! local extension ratio 
        !! osgd11(k,i) , osgd22(k,i), osdgd11dl(k,i), osdgd22dl(k,i)
        !! are defined and evaluated at rest state at inital 
        !! and stored in the module 
        js = dsqrt(osgu11(k,i) * osgu22(k,i) * gd11 * gd22)
        djsds = 1d0 / (2d0 * js) * ( (osdgu11ds(k,i) * osgu22(k,i) * gd11 * gd22) &
                                     + (osgu11(k,i) * osdgu22ds(k,i) * gd11 * gd22) &
                                     + (osgu11(k,i) * osgu22(k,i) * dgd11ds * gd22) &
                                     + (osgu11(k,i) * osgu22(k,i) * gd11 * dgd22ds) )                        
        
        Call SurfaceTensionDependOnCoverage(js, djsds, gu11, gu22, dgu22ds, b11, b22,&
                                           c112, c211, c222, normal_stress, tangent_stress)
!        Call SurfaceTensionConstant(gu11, gu22, dgu22ds, b11, b22,&
!                             c112, c211, c222, normal_stress, tangent_stress) 
!        Call ElasticStress(js, djsds, gu11, gu22, dgu22ds, osgu11(k,i), osgu22(k,i), &
!            osdgu22ds(k,i), b11, b22, c112, c211, c222, normal_stress, tangent_stress)                                    
                                           
                                           
        normal_stress1o_s(k,i) = normal_stress
        tangent_stress1o_s(k,i) = tangent_stress
                         
        dlds1o_s(k,i) = dlds       
        s1o_s(k,i) = s
        r1o_s(k,i) = rs
        z1o_s(k,i) = zs
        tr1o_s(k,i) = trs
        tz1o_s(k,i) = tzs
        d2rds2o1_s(k,i) = d2rds2
        d2zds2o1_s(k,i) = d2zds2
        gu22o1_s(k,i) = gu22
        dgu22dso1_s(k,i) = dgu22ds
     
        nr1o_s(k,i) = nrs
        nz1o_s(k,i) = nzs
        partofcovarintderiv1o_s(k,i) = dgu22ds + (c112+c222+c222) * gu22 + c211 * gu11
        meancurv1o_s(k,i) = b11 * gu11 + b22 * gu22
              
    ENDDO !! k-loop
    
ENDDO !! i-loop
        
END SUBROUTINE



SUBROUTINE for1_special_arrayo_cubs_gq()
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE

DOUBLE PRECISION :: nrs,nzs,trs,tzs,sfield,estress
DOUBLE PRECISION :: normal_force, tangent_force

INTEGER :: i,k,j    !!!  kth abscissa points,ith element;

Do i=1,N
 
    Do k = 1,2*sqm      
        
        sfield = 0d0
        Do j = 1,N+1
            sfield = car_cl1o_s(j,k,i) * dphidn1(j) + sfield
        ENDDO
        estress = sfield*sfield/2d0 * beta !beta -- electric capillary number
        
        !!! evaluate and store surface_force                
        nrs = nr1o_s(k,i) 
        nzs = nz1o_s(k,i)           
        trs = tr1o_s(k,i)
        tzs = tz1o_s(k,i)

!       normal_force = div_n - estress !+ b11 * t11 + b22 * t22 !including surface tension
!       normal_force = -(b11 * t11 + b22 * t22) ! neglect surface tension
        normal_force = - normal_stress1o_s(k,i) - estress !- div_n * kappa
        !write(11,*)  i, k, -(b11 * t11 + b22 * t22)
!        normal_force = - div_n - estress
        
        tangent_force = - tangent_stress1o_s(k,i)
        ! r component
        sforce_r1o_s(k,i) = nrs * normal_force + trs * tangent_force
        ! z component
        sforce_z1o_s(k,i) = nzs * normal_force + tzs * tangent_force
               
!       write(76,*) i,',',k,',',s,',',sforce_r1o(k,i),',',sforce_z1o(k,i)
    enddo !! for k 
enddo !! for i

END SUBROUTINE




SUBROUTINE geo1_metrics_arrayo_cubs_gq(istep)
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
INTEGER,INTENT(IN) :: istep
DOUBLE PRECISION :: s,rs,nrs,zs,nzs,dlds,trs,tzs
DOUBLE PRECISION :: temp,dzds,drds,d2zds2,d2rds2
DOUBLE PRECISION :: ps1i,ar11i,ar12i,ar13i,az11i,az12i,az13i,r1i,z1i
!! derivatives
!DOUBLE PRECISION :: drdl,dzdl,d2lds,d2rdl,d2zdl 
!! metric tensors for deformed state,g12 and g21 are zero
!! u means upper index, contravariant; d means lower index, covariant
DOUBLE PRECISION :: gu11,gu22,gd11,gd22
DOUBLE PRECISION :: dgu11ds,dgu22ds,dgd11ds,dgd22ds
!! christoffel symbol, c111, c122, c212 are zero
!! the first number is contravariant index and the last 2 numbers are covariant index
DOUBLE PRECISION :: c112,c211,c222
!! curvature tensor, b12 and b21 are zero
!! two numbers are covariant (lower) index
DOUBLE PRECISION :: b11,b22 
!! local extension ratio, js
DOUBLE PRECISION :: js, djsds
!! contravariant elastic stress resultants
!DOUBLE PRECISION :: t11, t22, dt22ds
!! surface tension, as a function of surface coverage
!DOUBLE PRECISION :: gamma
DOUBLE PRECISION :: normal_stress, tangent_stress

INTEGER :: i,k   !!!  kth abscissa points,ith element;

IF (mod(istep,10) == 1) THEN
   write(20,*) 'istep = ',istep, 's, zs, js(ratio of local area, deformed / undeformed)'
ENDIF

Do i=1,N
    ps1i = ps1(i) 
    ar11i = ar1(1,i)
    ar12i = ar1(2,i) 
    ar13i = ar1(3,i)
    az11i = az1(1,i)
    az12i = az1(2,i)
    az13i = az1(3,i)
    r1i = r1(i)
    z1i = z1(i)
    
    Do k=1,qm
!        IF (i == 39 .and. k == 3) THEN
!          write(*,*) i
!        ENDIF 
        s = s1o(k,i)
        temp = s-ps1i
        rs=ar11i*temp**3+ar12i*temp**2+ar13i*temp+r1i
        zs=az11i*temp**3+az12i*temp**2+az13i*temp+z1i
        !first order derivatives
        drds=3d0*ar11i*temp**2+2d0*ar12i*temp+ar13i
        dzds=3d0*az11i*temp**2+2d0*az12i*temp+az13i
        dlds = dsqrt(drds**2+dzds**2)
!            !also unit tangent vector r,z component
!        drdl = drds / dlds 
!        dzdl = dzds / dlds
       
        !second order derivatives
        d2rds2 = 6d0*ar11i*temp+2d0*ar12i
        d2zds2 = 6d0*az11i*temp+2d0*az12i
             
!        d2lds = (drds * d2rds + dzds * d2zds) / dlds        
!        d2rdl = d2rds / (dlds*dlds) - drds * d2lds / (dlds**3)
!        d2zdl = d2zds / (dlds*dlds) - dzds * d2lds / (dlds**3)
 
        !normal vector r,z component,pointing outward
        nrs= dzds/dlds
        nzs= -drds/dlds
        ! a_2 vector r,z component        
        trs = drds
        tzs = dzds
        
        !For axisymetric curvature
!        drdz = drds/dzds
!         !second derivative of r in z
!        d2rdz = d2rds/(dzds*dzds)-d2zds*drds/(dzds**3)        
!        temp1 = dsqrt(1d0+drdz*drdz)
!        ! div of unit normal vector
!        ! equal to twice the mean curvature
!        div_n = 1d0/(rs*temp1)-1d0/temp1**3*d2rdz
!        div_n = sign(div_n,drds*d2zds-dzds*d2rds)
        
        ! metric tensors
        gu11 = 1d0 / (rs*rs)
        gd11 = rs * rs
        gu22 = 1d0 / (drds**2 + dzds**2)        
        gd22 = drds**2 + dzds**2
        
        !dervivative of metric tensors
        dgd11ds = 2d0 * rs * drds
        dgd22ds = 2d0 * (drds * d2rds2 + dzds * d2zds2)
        dgu11ds = -2d0 * drds / rs**3
        dgu22ds = -2d0 * (drds * d2rds2 + dzds * d2zds2) / (drds**2 + dzds**2)**2 
        
        
        ! Christoffel symbols
        c112 = drds / rs
        c211 = -rs * drds / (drds**2 + dzds**2)
        c222 = (drds * d2rds2 + dzds * d2zds2) / (dzds**2 + drds**2)
        
        ! Curvature tensors
        b11 = - dzds * rs / dsqrt(dzds**2 + drds**2) 
        b22 = (dzds * d2rds2 - drds * d2zds2) / dsqrt(dzds**2 + drds**2)

        !! local extension ratio 
        !! ogd11(k,i) , ogd22(k,i), odgd11dl(k,i), odgd22dl(k,i)
        !! are defined and evaluated at rest state at inital 
        !! and stored in the module 
        js = dsqrt(ogu11(k,i) * ogu22(k,i) * gd11 * gd22)
        djsds = 1d0 / (2d0 * js) * ( (odgu11ds(k,i) * ogu22(k,i) * gd11 * gd22) &
                                     + (ogu11(k,i) * odgu22ds(k,i) * gd11 * gd22) &
                                     + (ogu11(k,i) * ogu22(k,i) * dgd11ds * gd22) &
                                     + (ogu11(k,i) * ogu22(k,i) * gd11 * dgd22ds) )        
        
        Call SurfaceTensionDependOnCoverage(js, djsds, gu11, gu22, dgu22ds, b11, b22,&
                                            c112, c211, c222, normal_stress, tangent_stress) 
!        Call SurfaceTensionConstant(gu11, gu22, dgu22ds, b11, b22,&
!                             c112, c211, c222, normal_stress, tangent_stress) 
!        Call ElasticStress(js, djsds, gu11, gu22, dgu22ds, ogu11(k,i), ogu22(k,i), &
!            odgu22ds(k,i), b11, b22, c112, c211, c222, normal_stress, tangent_stress) 
                                            
        normal_stress1o(k,i) = normal_stress
        tangent_stress1o(k,i) = tangent_stress
        
        dlds1o(k,i) = dlds
        s = s1o(k,i)
        r1o(k,i) = rs
        z1o(k,i) = zs
        tr1o(k,i) = trs
        tz1o(k,i) = tzs
        nr1o(k,i) = nrs
        nz1o(k,i) = nzs
        d2rds2o1(k,i) = d2rds2
        d2zds2o1(k,i) = d2zds2
        gu11o1(k,i) = gu11
        gu22o1(k,i) = gu22
        dgu22dso1(k,i) = dgu22ds
        
        partofcovarintderiv1o(k,i) = dgu22ds + (c112+c222+c222) * gu22 + c211 * gu11
        meancurv1o(k,i) = b11 * gu11 + b22 * gu22
       
        IF (mod(istep,10) == 1) THEN
             write(20,400) s,',',zs,',',js
        ENDIF
    enddo !! for k 
enddo !! for i

400 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)
END SUBROUTINE


SUBROUTINE for1_arrayo_cubs_gq()
use mod_elasbody
use mod_gaussqxw
IMPLICIT NONE

DOUBLE PRECISION :: nrs,nzs,trs,tzs,sfield,estress
DOUBLE PRECISION :: normal_force, tangent_force

INTEGER :: i,k,j    !!!  kth abscissa points,ith element;

Do i=1,N

    Do k=1,qm
       
        sfield = 0d0
        Do j = 1,N+1
            sfield = car_cl1o(j,k,i) * dphidn1(j) + sfield
        Enddo
        estress = sfield*sfield/2d0 * beta !beta -- electric capillary number
        !!! evaluate and store surface_force 
               
        nrs = nr1o(k,i) 
        nzs = nz1o(k,i)           
        trs = tr1o(k,i)
        tzs = tz1o(k,i)

!        normal_force = div_n - estress !+ b11 * t11 + b22 * t22 !including surface tension
!       normal_force = -(b11 * t11 + b22 * t22) ! neglect surface tension
        normal_force = - normal_stress1o(k,i) - estress !- div_n * kappa
        !write(11,*)  i, k, -(b11 * t11 + b22 * t22)
!        normal_force = - div_n - estress
        
        tangent_force = - tangent_stress1o(k,i)
        ! r component
        sforce_r1o(k,i) = nrs * normal_force + trs * tangent_force
        ! z component
        sforce_z1o(k,i) = nzs * normal_force + tzs * tangent_force
               
!       write(76,*) i,',',k,',',s,',',sforce_r1o(k,i),',',sforce_z1o(k,i)
    enddo !! for k 
enddo !! for i

END SUBROUTINE
